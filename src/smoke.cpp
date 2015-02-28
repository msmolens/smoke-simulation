// Max Smolens
// max@cs.unc.edu
//
// COMP 259 Final Project
// May 2003
//
// Smoke Simulation

#include <particle/papi.h>
#include <particle/p_vector.h>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <list>
#include <vector>
#include <GL/glut.h>
#include <glui.h>
#include "ImageWriter.h"
#include "turbulence.h"
#include "cdist.h"

using namespace std;

#ifdef WIN32
#pragma warning (disable:4305) /* disable bogus conversion warnings */
#define drand48() (((float) rand())/((float) RAND_MAX))
#define lrand48() ((rand() << 16) ^ rand())
#define srand48(x) srand(x)
#endif

// #define DEPTH_TEST

#define GL_ASSERT() {GLenum sci_err; while ((sci_err = glGetError()) != GL_NO_ERROR) \
			cerr << "OpenGL error: " << (char *)gluErrorString(sci_err) << " at " << __FILE__ <<":" << __LINE__ << endl;}

// gui
GLUI *glui = NULL;
int main_window;
enum glui_controls {
	DRAW_GROUND_CHECKBOX,
	DRAW_AXES_CHECKBOX,
	DRAW_GLOBAL_FORCE_CHECKBOX,
	DO_MOTION_CHECKBOX,
	PAUSE_CHECKBOX,
	SCREENDUMP_BUTTON,
	RENDER_METHOD_RADIOGROUP,
	QUAD_SIZE_SPINNER,
	PARTICLE_SPREAD_SPINNER,
	CCONV1_EDITTEXT,
	CCONV2_EDITTEXT,
	TSIM_EDITTEXT,
	CCOOL_EDITTEXT,
	CFRIC1_EDITTEXT,
	CFRIC2_EDITTEXT,
	CDISS_EDITTEXT,
	NUM_NEW_PARTICLES_SPINNER,
	SMOKE_COLOR_SPINNER,
	NEW_CLUSTER_DELAY_SPINNER,
	GLOBAL_FORCE_RESET_BUTTON,
	GLOBAL_FORCE_START_BUTTON,
	GLOBAL_FORCE_STRENGTH1_EDITTEXT,
	GLOBAL_FORCE_STRENGTH2_EDITTEXT,
	GLOBAL_FORCE_DIRECTION1_XY_TRANSLATION,
	GLOBAL_FORCE_DIRECTION1_Z_TRANSLATION,
	GLOBAL_FORCE_DIRECTION2_XY_TRANSLATION,
	GLOBAL_FORCE_DIRECTION2_Z_TRANSLATION,
	GLOBAL_FORCE_DURATION_EDITTEXT

};
GLUI_EditText *cconv1_edittext, *cconv2_edittext, *tsim_edittext, *ccool_edittext, *cfric1_edittext, *cfric2_edittext, *cdiss_edittext, *global_force_strength1_edittext, *global_force_strength2_edittext, *global_force_duration_edittext;
GLUI_Checkbox *draw_ground_checkbox, *draw_axes_checkbox, *do_motion_checkbox, *pause_checkbox, *draw_global_force_checkbox;
GLUI_RadioGroup *render_method_radiogroup;
GLUI_Button *screendump_button, *global_force_start_button, *global_force_reset_button;
GLUI_Spinner *quad_size_spinner, *particle_spread_spinner, *num_new_particles_spinner, *smoke_color_spinner, *new_cluster_delay_spinner;
GLUI_StaticText *cluster_statictext;
GLUI_Translation *global_force_direction1_xy_translation, *global_force_direction2_xy_translation, *global_force_direction1_z_translation, *global_force_direction2_z_translation;;

// booleans
bool double_buffer = true, antialias = true;
bool fullscreen = false;
int draw_axes = 0, draw_ground = 0, draw_global_force = 0, freeze = 0, do_motion = 0;
int num_steps = 1, spot_tex_id = -1;
int render_method = 0;

// window size
static int w = 512;
static int h = 512;

// simulation constants
static float cconv1       = 3000;
static float cconv2       = 3000;
static const float t0sim  = 200;
static float tsim         = 70;
static float ccool        = 500;
static float cfric1       = 1.0;
static float cfric2       = 1.5;
static float cdiss        = 0.10;

// global force
bool global_force_applied = false;
static int global_force_start_ticks = 0;
static int global_force_stop_ticks  = 0;
static int global_force_duration    = 500;
static const float GLOBAL_FORCE_STRENGTH_SCALE = 0.0001;
static float global_force_strength1  = 2.0;
static float global_force_strength2  = 1.0;
static pVector global_force_direction1(1.0, 1.0, 0.0);
static pVector global_force_direction2(-1.0, -1.0, 0.0);
static pVector global_force1(0.0, 0.0, 0.0);
static pVector global_force2(0.0, 0.0, 0.0);
static pVector gf_draw(0.0, 0.0, 0.0);
static pVector gf(0.0, 0.0, 0.0);

// rendering constants
static const int MAX_PARTICLES         = 100;
static float particle_spread           = 0.50;
static const float PARTICLE_SPREAD_MIN = 0.15;
static const float PARTICLE_SPREAD_MAX = 1.25;
static float quad_size                 = 0.5;
static const float QUAD_SIZE_MIN       = 0.01;
static int num_new_particles           = 10;
static const float LOWEST_DISPLAYED_DENSITY = 0.05;
static float smoke_color               = 0.65;
static int new_cluster_delay          = 1;

// particle groups list
list<int> groups;

// simulation clocks
static uint ticks       = 0;
static uint movie_ticks = 0;

// viewing parameters
static pVector eye(0.0, -10.0, 5.0);
static pVector center(0.5, 0.5, 5.0);
static const pVector up(0.0, 0.0, 1.0);
static float rot       = 0.0f;
static float view_dist = 10.0f;

// source position list (turbulence space)
static const int NUM_SOURCES = 5;
int cur_source = 0;
float ss = 30.0;			// source position scale
float sources[] = {
	0.0, 0.0, 0.0 ,
	1.0*ss, 1.0*ss, 0.0,
	1.0*ss, -1.0*ss, 0.0,
	-1.0*ss, 1.0*ss, 0.0,
	-1.0*ss, -1.0*ss, 0.0
};

// image writer
static ImageWriter writer;
static bool writing = false;

// glut font
static const int font=(int)GLUT_BITMAP_8_BY_13;


//////// BEGIN FUNCTIONS ////////

// render a string
void
render_string(float x, float y, void *font, char *string)
{
	char *c;
	glRasterPos2f(x, y);
	for (c=string; *c != '\0'; c++) {
		glutBitmapCharacter(font, *c);
	}
}

// set random values
void
set_rand3_zpos(float *r, float scale = 1.0)
{
	r[0] = drand48()*scale*(drand48()<0.5 ? 1.0 : -1.0);
	r[1] = drand48()*scale*(drand48()<0.5 ? 1.0 : -1.0);
	r[2] = drand48()*scale;	// z will be positive
}
void
set_rand3(float *r, float scale = 1.0)
{
	r[0] = drand48()*scale*(drand48()<0.5 ? 1.0 : -1.0);
	r[1] = drand48()*scale*(drand48()<0.5 ? 1.0 : -1.0);
	r[2] = drand48()*scale*(drand48()<0.5 ? 1.0 : -1.0);
}

// symmetric gaussian centered at origin
// from PSpray demo
inline float
Gaussian2(float x, float y, float sigma)
{
// the sqrt of 2 pi.
#define SQRT2PI 2.506628274631000502415765284811045253006
	return exp(-0.5 * (x*x + y*y) / (sigma*sigma)) / (SQRT2PI * sigma);
}

// adapted from PSpray demo
void
MakeGaussianSpotTexture()
{
	static const int DIM = 32;
	static const int DIM2 = (DIM>>1);

	glGenTextures(1, (GLuint *)&spot_tex_id);
	glBindTexture(GL_TEXTURE_2D, spot_tex_id);

	float *img = new float[DIM*DIM];

	for(int y=0; y<DIM; y++)
	{
		for(int x=0; x<DIM; x++)
		{
			if(x==0 || x==DIM-1 || y==0 || y==DIM-1)
				img[y*DIM+x] = 0;
			else
			{
//				img[y*DIM+x] = 5.0 * Gaussian2(x-DIM2, y-DIM2, (DIM*0.2));
				img[y*DIM+x] = 1.0 * Gaussian2(x-DIM2, y-DIM2, (DIM*0.2));
			}
		}
	}
	
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	float col[4] = {1.f, 1.f, 1.f, 1.f};
	glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, col);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	GL_ASSERT();

	gluBuild2DMipmaps(GL_TEXTURE_2D, GL_ALPHA8, DIM, DIM, GL_ALPHA, GL_FLOAT, img);

	delete [] img;

}


// draw as a splat texture on a quad
// adapted from PSpray demo
void
draw_textured_quad(const pVector &view, const pVector &up, float size_scale = 1.0f)
{
	int cnt = pGetGroupCount();
	if (cnt < 1)
		return;
	
	pVector *ppos = new pVector[cnt];
	pVector *size = new pVector[cnt];
	pGetParticles(0, cnt, (float *)ppos, NULL, NULL, (float *)size);

	// compute the vectors from the particle to the corners of its quad
	//         the particle is at the center of the x
	// 3-2     V0, V1, V2 and V3 go from there to the vertices
	// |x|     the texcoords are (0,0), (1,0), (1,1), and (0,1) respectively
	// 0-1     we clamp the texture so the rest is transparent

	pVector right = view ^ up;
	right.normalize();
	pVector nup = right ^ view;
	right *= size_scale;
	nup   *= size_scale;

	pVector V0 = -(right + nup);
	pVector V1 = right - nup;
	pVector V2 = right + nup;
	pVector V3 = nup - right;

	glBegin(GL_QUADS);

	for (int i = 0; i < cnt; i++) {
		// don't render if density is low
		// particles should really already have been removed from the system
		// but particle system api doesn't have that type of control
		// this works well enough, so I didn't add that to the api
		if (size[i].x < LOWEST_DISPLAYED_DENSITY)
			continue;

 		pVector &p = ppos[i];
		glColor4f(smoke_color, smoke_color, smoke_color, size[i].x);
		
 		glTexCoord2f(0,0);
 		pVector ver = p + V0;
 		glVertex3fv((GLfloat *)&ver);

 		glTexCoord2f(1,0);
 		ver = p + V1;
 		glVertex3fv((GLfloat *)&ver);

 		glTexCoord2f(1,1);
 		ver = p + V2;
 		glVertex3fv((GLfloat *)&ver);

 		glTexCoord2f(0,1);
 		ver = p + V3;
 		glVertex3fv((GLfloat *)&ver);
	}
	
	glEnd();

	delete [] ppos;
	delete [] size;
}


// draws groups in sorted back to front order
void
draw_textured_groups()
{
	vector<cdist> gdraw;
	pVector ppos;

	// calculate distance to each group
	for (list<int>::iterator it = groups.begin(); it != groups.end(); ++it) { 
		pCurrentGroup((*it));
		if (pGetGroupCount() < 1)
			continue;
		pGetParticles(0, 1, (float *)&ppos);
		pVector v = ppos - eye;
		float dist = v.length();
		gdraw.push_back(cdist((*it), dist));
	}

	// sort by distance to camera
	sort(gdraw.begin(), gdraw.end());

	// update view vector
	pVector view = center - eye;
	view.normalize();

	glColor4f(0.0,0.2,0.3,1.0);

	glEnable(GL_TEXTURE_2D);
	// draw groups in back to front order
	vector<cdist>::reverse_iterator it;
	for (it = gdraw.rbegin(); it != gdraw.rend();  ++it) {
		pCurrentGroup(it->g);
		draw_textured_quad(view, up, quad_size);
	}

	glDisable(GL_TEXTURE_2D);
}


// T(t) as described in Holtkamper paper
// similar to Newtonians law of cooling
float
cooling(float t)
{
	return tsim + (t0sim - tsim) * exp(-t/ccool);
}

// adds new particle cluster to simulation
// s: source index
// g: group number (don't specify to create new group)
// returns new group id
int
new_group(int s, int g = 0)
{
	float o[3], v[3];
	float *p = &sources[3*s];
	static const float vel_scale = 0.001;

	// create new particle group and set active
	if (g == 0)
		g = pGenParticleGroups(1, MAX_PARTICLES);
	pCurrentGroup(g);

	// set random velocity (uniform for all particles)
	set_rand3_zpos(v, vel_scale);
	pVelocity(v[0], v[1], v[2]);

	// move source position in turbulence space opposite velocity
	p[0] -= v[0]*1000.0;
	p[1] -= v[1]*1000.0;
	p[2] -= v[2]*1000.0;

	float fp[3];
	fp[0] = (p[0] - floor(p[0]))*1.5;
	fp[1] = (p[1] - floor(p[1]))*1.5;
	fp[2] = (p[2] - floor(p[2]))*1.5;

	// add particles to current group
	for (int i = 0 ; i < num_new_particles; i ++) {

		// offset (in both turbulence space and simulation space)
		set_rand3(o, particle_spread);
		o[2] *= 0.5;

		// density
		float d = turbulence(p[0]+o[0], p[1]+o[1], p[2]+o[2]);
		pSize(d);

		// color for point rendering
		pColor(1.0, 1.0, 1.0, d);

		// add particle
		pVertex(fp[0]+o[0], fp[1]+o[1], fp[2]+o[2]);
	}

	return g;
}

// runs a timestep of smoke simulation
void
step_smoke()
{
	// pick random source
//	cur_source = (cur_source + 1) % NUM_SOURCES;
	if (ticks % new_cluster_delay == 0) {
		cur_source = (int)((float)NUM_SOURCES * rand()/(RAND_MAX+1.0));
		groups.push_back(new_group(cur_source));
	}

	// check if global force applied
	if (global_force_applied) {
		// calculate smooth curve between initial force and ending force
		// see paper for formula
		float a = float(ticks - global_force_start_ticks)/float(global_force_duration);
		a = a * a;
		gf_draw = (global_force2 - global_force1) * (3.0 * a - 2.0 * a * a) + global_force1;
		gf = gf_draw * GLOBAL_FORCE_STRENGTH_SCALE;
		if (ticks >= global_force_stop_ticks) {
			global_force_applied = false;
			gf      = pVector(0.0, 0.0, 0.0);
			global_force_start_button->enable();
			global_force_reset_button->enable();
		}
	}
	
	// update groups
	list<int>::iterator it;
	for (it = groups.begin(); it != groups.end(); ++it) { 

		pCurrentGroup((*it));

		// delete empty groups
 		int cnt = pGetGroupCount();
 		if (cnt < 1) {
			pDeleteParticleGroups(*it);
			groups.erase(it++);
			continue;
		}

		// update particle parameters
		pVector *psize = new pVector[cnt]; // size
		pVector *pvel  = new pVector[cnt]; // velocities
		float *pages   = new float[cnt*3]; // ages
		
 		// args: pos color vel size age
  		pGetParticles(0, cnt, NULL, NULL, (float *)pvel, (float *)psize, pages);

		// calculate new velocity based on forcse
		pVector *ovel = &pvel[0]; // original velocity
		pVector vel   = pvel[0];  // new velocity
		
		// gravitational force
		pVector g(0.0, 0.0, -0.00018);
		
		// convection force
		float conv = (cconv1/tsim) - (cconv2/cooling(pages[0]));
		conv *= 0.00001;
		pVector c(0.0, 0.0, conv);

		// frictional force
		float normv = ovel->length();
		float fric  = -cfric1 * pow(normv, cfric2);
		pVector f = *ovel;
		f.normalize();
		f *= fric * 0.0005;

		// sum all forces
		vel += g + c + f + gf;
		pRandomVelocity(PDPoint, vel.x, vel.y, vel.z);
		
		// densities (stored in size.x)
		for (int i = 0; i < cnt; i++) {
			float d = psize[i].x - (cdiss * normv);
			d -= psize[i].x / 10000.0; // magic number
		                        // should really be remaining lifetime of particle
                                        // (but I don't constrain age)
			if (d < 0.0)
				d = 0.0;
			psize[i].x = d;
		}

		pSetSizes((float *)psize);

		// update positions
		pMove();

                // kill particles that start to fall
                pSinkVelocity(false, PDBox,
                              -P_MAXFLOAT, -P_MAXFLOAT, -0.01,
                              P_MAXFLOAT, P_MAXFLOAT, P_MAXFLOAT);

		delete [] psize;
		delete [] pvel;
		delete [] pages;
	}

	ticks++;
}

// main drawing function
void
draw()
{
	if (glutGetWindow() != main_window) 
		glutSetWindow(main_window);  

	// update cluster count
	char c[25];
	snprintf(c, 25, "Clusters: %d", groups.size());
	cluster_statictext->set_text(c);
	
	glLoadIdentity();
	
	if (spot_tex_id < 0)
		MakeGaussianSpotTexture();

#ifdef DEPTH_TEST
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
#else
	glClear(GL_COLOR_BUFFER_BIT);
#endif

	// auto-rotation
	if (do_motion) {
		rot += 0.01;
		if (rot >= 2*M_PI)
			rot -= 2*M_PI;
	}
	
	// calculate eye position
	eye.x = view_dist * sin(rot);
	eye.y = view_dist * cos(rot);

	gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z);

	// draw ground
	glColor4f(0.0,0.3,0.4,1.0);
	if (draw_ground) {
		glBegin(GL_QUADS);
		glVertex3f(-2,-2,0);
		glVertex3f(-2,3,0);
		glColor4f(0.0,0.32,0.4,1.0);
		glVertex3f(3,3,0);
		glVertex3f(3,-2,0);
		glEnd();
	}

	// step simulation
	if (!freeze) {
		for(int step = 0; step < num_steps; step++) {
			step_smoke();
		}
	}
	
	// render smoke
	if (render_method == 1) { // draw as points
		list<int>::iterator it;
		for (it = groups.begin(); it != groups.end(); ++it) { 
			pCurrentGroup((*it));
			pDrawGroupp(GL_POINTS);
		}
	} else {		// draw as textured quads
		draw_textured_groups();
	}

	// draw axes
	if (draw_axes) {
		glPushMatrix();
		glTranslatef(center.x, center.y, 4.0);
		glBegin(GL_LINES);
		glColor3f(1.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(1.0, 0.0, 0.0);
		glColor3f(0.0, 1.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 1.0, 0.0);
		glColor3f(0.0, 0.0, 1.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, 1.0);
		glEnd();
		glPopMatrix();
	}

	// draw global force
	if (draw_global_force) {
		if (global_force_applied) {
			glPushMatrix();
			glTranslatef(center.x, center.y, 4.0);
			glBegin(GL_LINES);
			glColor3f(1.0, 1.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3fv(&gf_draw.x);
			glEnd();
			glPopMatrix();
		} else {
			glPushMatrix();
			glTranslatef(center.x, center.y, 4.0);
			glBegin(GL_LINES);
			glColor3f(0.0, 1.0, 1.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3fv(&global_force_direction1.x);
			glColor3f(1.0, 0.0, 1.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3fv(&global_force_direction2.x);
			glEnd();
			glPopMatrix();
		}
	}

	// render text
#if 0
	glPushMatrix();
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        gluOrtho2D(0, w, 0, h);
        glScalef(1.0, -1.0, 1.0);
        glTranslatef(0.0, -h, 0.0);
        glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	char str[50];
	snprintf(str, 50, "Clusters: %d", groups.size());
	render_string(5, 20, (void *)font, str);
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
#endif

	GL_ASSERT();

	if (double_buffer)
		glutSwapBuffers();
	glutPostRedisplay();

	// dump images
	if (writing) {
		writer.output(movie_ticks);
		movie_ticks++;
	}
}

// reshape callback
void reshape(int _w, int _h)
{
	w = _w;
	h = _h;

	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(40, w / double(h), 1, 100);
	glMatrixMode(GL_MODELVIEW);
	
#ifdef DEPTH_TEST
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
#else
	glClear(GL_COLOR_BUFFER_BIT);
#endif

	writer.set_dimensions(w, h);
}

// menu handler
void
menu(int item)
{
	static int OldWidth, OldHeight;

	switch(item)
	{
	case 'a':
		rot += 0.05;
		if (rot >= 2*M_PI)
			rot -= 2*M_PI;
		break;
	case 'd':
		rot -= 0.05;
		if (rot <= 0.0)
			rot += 2*M_PI;
		break;
	case 'w':
		view_dist -= 0.1;
		break;
	case 's':
		view_dist += 0.1;
		break;
	case 't':
		center.z += 0.05;
		eye.z += 0.05;
		break;
	case 'g':
		center.z -= 0.05;
		eye.z -= 0.05;
		if (center.z < 2.0 || eye.z < 2.0) {
			center.z = 2.0;
			eye.z = 2.0;
		}
		break;
	case '/':
		particle_spread -= 0.01;
		if (particle_spread < PARTICLE_SPREAD_MIN)
			particle_spread = PARTICLE_SPREAD_MIN;
		break;
	case '*':
		particle_spread += 0.01;
		if (particle_spread > PARTICLE_SPREAD_MAX)
			particle_spread = PARTICLE_SPREAD_MAX;
		break;
	case 'i':
		writing = !writing;
		cout << "writing: " << writing << endl;
		break;
	case 'r':
		antialias = !antialias;
		if(antialias)
		{
			glEnable(GL_LINE_SMOOTH);
			glEnable(GL_POINT_SMOOTH);
		}
		else
		{
			glDisable(GL_LINE_SMOOTH);
			glDisable(GL_POINT_SMOOTH);
		}
		break;
	case 'c':
		do_motion = !do_motion;
		break;
	case 'f':
		fullscreen = !fullscreen;
		if (fullscreen) {
			OldWidth = glutGet(GLenum(GLUT_WINDOW_WIDTH));
			OldHeight = glutGet(GLenum(GLUT_WINDOW_HEIGHT));
			glutSetCursor(GLUT_CURSOR_NONE);
			glutFullScreen(); 
		} else	{
			glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
			glutReshapeWindow(OldWidth, OldHeight);
		}
		break;
	case 'p':
		render_method = 0;
		glDisable(GL_TEXTURE_2D);
		break;
	case 'o':
		render_method = 1;
		glEnable(GL_TEXTURE_2D);
		break;
	case 'h':
		draw_ground = !draw_ground;
		break;
	case 'k':
		draw_axes = !draw_axes;
		break;
	case 'x':
		freeze = !freeze;
		break;
	case '+':
		quad_size += 0.01;
		break;
	case '-':
		quad_size -= 0.01;
		if (quad_size < QUAD_SIZE_MIN) {
			quad_size = QUAD_SIZE_MIN;
		}
		break;
	case 'q':
	case 'Q':
	case 27:
		exit(0);
		break;
	}
	
	if (item > '0' && item <= '9') {
		num_steps = item - '0';
		pTimeStep(1 / float(num_steps));
	}
	
	glui->sync_live();

	glutPostRedisplay();
}

// keypress callback
void
key_press(unsigned char key, int x, int y)
{
	menu((int) key);
}

// glui callback
void
glui_cb(int id)
{
	switch(id) {
	case SCREENDUMP_BUTTON:
		writing = !writing;
		cout << "writing: " << writing << endl;
		break;
	case RENDER_METHOD_RADIOGROUP:
		render_method = (render_method > 0) ? 1 : 0;
		break;
	case GLOBAL_FORCE_DIRECTION1_XY_TRANSLATION:
	case GLOBAL_FORCE_DIRECTION1_Z_TRANSLATION:
		global_force_direction1.normalize();
		break;
	case GLOBAL_FORCE_DIRECTION2_XY_TRANSLATION:
	case GLOBAL_FORCE_DIRECTION2_Z_TRANSLATION:
		global_force_direction2.normalize();
		break;
	case GLOBAL_FORCE_START_BUTTON:
		global_force_start_button->disable();
		global_force_reset_button->disable();
		global_force_start_ticks = ticks;
		global_force_stop_ticks  = ticks + global_force_duration;
		global_force1 = global_force_direction1 * global_force_strength1;
		global_force2 = global_force_direction2 * global_force_strength2; 
		global_force_applied = true;
		break;
	case GLOBAL_FORCE_RESET_BUTTON:
		global_force_direction1 = pVector(1.0, 1.0, 0.0);
		global_force_direction2 = pVector(-1.0, -1.0, 0.0);
		global_force_direction1.normalize();
		global_force_direction2.normalize();
		break;
	default:
		break;
	}

	glui->sync_live();
}

// glui initialization
void
init_glui()
{
	glui = GLUI_Master.create_glui("Controls", 0, 50, 50);
	glui->set_main_gfx_window(main_window);
	GLUI_Master.set_glutIdleFunc(draw);
	GLUI_Master.set_glutReshapeFunc(reshape);

	GLUI_Panel *panel0 = glui->add_panel("", GLUI_PANEL_NONE);
	panel0->set_alignment(GLUI_ALIGN_LEFT);

	GLUI_Panel *panel1 = glui->add_panel_to_panel(panel0, "Simulation constants");
	cconv1_edittext =
		glui->add_edittext_to_panel(panel1, "cconv1", GLUI_EDITTEXT_FLOAT, &cconv1,
				   CCONV1_EDITTEXT);
	cconv1_edittext->set_float_limits(2500.0, 3500.0);
 	cconv2_edittext =
		glui->add_edittext_to_panel(panel1, "cconv2", GLUI_EDITTEXT_FLOAT, &cconv2,
				   CCONV2_EDITTEXT, glui_cb);
	cconv2_edittext->set_float_limits(2500.0, 3500.0);
 	tsim_edittext =
		glui->add_edittext_to_panel(panel1, "tsim", GLUI_EDITTEXT_FLOAT, &tsim,
				   TSIM_EDITTEXT);
	tsim_edittext->set_float_limits(0.0, 200.0);
 	ccool_edittext =
		glui->add_edittext_to_panel(panel1, "ccool", GLUI_EDITTEXT_FLOAT, &ccool,
				   CCOOL_EDITTEXT);
	ccool_edittext->set_float_limits(1.0, 1000.0);
 	cfric1_edittext =
		glui->add_edittext_to_panel(panel1, "cfric1", GLUI_EDITTEXT_FLOAT, &cfric1,
				   CFRIC1_EDITTEXT);
	cfric1_edittext->set_float_limits(0.1, 1.5);
 	cfric2_edittext =
		glui->add_edittext_to_panel(panel1, "cfric2", GLUI_EDITTEXT_FLOAT, &cfric2,
				   CFRIC2_EDITTEXT);
	cfric2_edittext->set_float_limits(0.1, 1.5);
 	cdiss_edittext =
		glui->add_edittext_to_panel(panel1, "cdiss", GLUI_EDITTEXT_FLOAT, &cdiss,
				   CDISS_EDITTEXT);
	cdiss_edittext->set_float_limits(0.0001, 0.5);
	glui->add_column_to_panel(panel0, false);

	GLUI_Panel *panel2 = glui->add_panel_to_panel(panel0, "Rendering");
	render_method_radiogroup =
		glui->add_radiogroup_to_panel(panel2, &render_method,
				     RENDER_METHOD_RADIOGROUP, glui_cb);
	glui->add_radiobutton_to_group(render_method_radiogroup, "Textured");
	glui->add_radiobutton_to_group(render_method_radiogroup, "Points");
	glui->add_separator_to_panel(panel2);
	draw_ground_checkbox =
		glui->add_checkbox_to_panel(panel2, "Draw ground", &draw_ground, DRAW_GROUND_CHECKBOX);
	draw_ground_checkbox->set_alignment(GLUI_ALIGN_LEFT);
	draw_axes_checkbox =
		glui->add_checkbox_to_panel(panel2, "Draw axes", &draw_axes, DRAW_AXES_CHECKBOX);
	draw_axes_checkbox->set_alignment(GLUI_ALIGN_LEFT);
	draw_global_force_checkbox =
		glui->add_checkbox_to_panel(panel2, "Draw global force", &draw_global_force, DRAW_GLOBAL_FORCE_CHECKBOX);
	draw_global_force_checkbox->set_alignment(GLUI_ALIGN_LEFT);
	do_motion_checkbox =
		glui->add_checkbox_to_panel(panel2, "Rotate", &do_motion, DO_MOTION_CHECKBOX);
	do_motion_checkbox->set_alignment(GLUI_ALIGN_LEFT);
	pause_checkbox =
		glui->add_checkbox_to_panel(panel2, "Pause", &freeze, PAUSE_CHECKBOX);
	pause_checkbox->set_alignment(GLUI_ALIGN_LEFT);

	glui->add_column_to_panel(panel2, false);
	panel2->set_alignment(GLUI_ALIGN_LEFT);
	quad_size_spinner =
		glui->add_spinner_to_panel(panel2, "Texture size", GLUI_SPINNER_FLOAT, &quad_size,
				  QUAD_SIZE_SPINNER);
	quad_size_spinner->set_float_limits(QUAD_SIZE_MIN, 3.0);
	quad_size_spinner->set_alignment(GLUI_ALIGN_LEFT);
	smoke_color_spinner =
		glui->add_spinner_to_panel(panel2, "Smoke color", GLUI_SPINNER_FLOAT, &smoke_color, SMOKE_COLOR_SPINNER);
	smoke_color_spinner->set_float_limits(0.0, 1.0);
	smoke_color_spinner->set_alignment(GLUI_ALIGN_LEFT);
	particle_spread_spinner =
		glui->add_spinner_to_panel(panel2, "Particle spread", GLUI_SPINNER_FLOAT, &particle_spread,PARTICLE_SPREAD_SPINNER);
	particle_spread_spinner->set_float_limits(0.15, 1.25);
	particle_spread_spinner->set_alignment(GLUI_ALIGN_LEFT);
 	num_new_particles_spinner =
		glui->add_spinner_to_panel(panel2, "Particles/cluster", GLUI_SPINNER_INT, &num_new_particles, NUM_NEW_PARTICLES_SPINNER);
	num_new_particles_spinner->set_int_limits(0, MAX_PARTICLES);
	num_new_particles_spinner->set_alignment(GLUI_ALIGN_LEFT);
 	new_cluster_delay_spinner =
		glui->add_spinner_to_panel(panel2, "New cluster delay", GLUI_SPINNER_INT, &new_cluster_delay, NEW_CLUSTER_DELAY_SPINNER);
	new_cluster_delay_spinner->set_int_limits(1, 15);
	new_cluster_delay_spinner->set_alignment(GLUI_ALIGN_LEFT);

	GLUI_Panel *panel6 = glui->add_panel("", GLUI_PANEL_NONE);

	GLUI_Panel *panel5 = glui->add_panel_to_panel(panel6, "Global force");
	panel5->set_alignment(GLUI_ALIGN_LEFT);
 	global_force_strength1_edittext =
		glui->add_edittext_to_panel(panel5, "Strength 1", GLUI_EDITTEXT_FLOAT, &global_force_strength1, GLOBAL_FORCE_STRENGTH1_EDITTEXT);
	global_force_strength1_edittext->set_float_limits(1.0, 10.0);
 	global_force_strength2_edittext =
		glui->add_edittext_to_panel(panel5, "Strength 2", GLUI_EDITTEXT_FLOAT, &global_force_strength2,GLOBAL_FORCE_STRENGTH2_EDITTEXT);
	global_force_strength2_edittext->set_float_limits(1.0, 10.0);
 	global_force_duration_edittext =
		glui->add_edittext_to_panel(panel5, "Duration", GLUI_EDITTEXT_INT, &global_force_duration, GLOBAL_FORCE_DURATION_EDITTEXT);
	global_force_duration_edittext->set_int_limits(1, 1000);
	global_force_start_button =
		glui->add_button_to_panel(panel5, "Apply global force", GLOBAL_FORCE_START_BUTTON, glui_cb);
	global_force_start_button->set_alignment(GLUI_ALIGN_CENTER);
	global_force_reset_button =
		glui->add_button_to_panel(panel5, "Reset directions", GLOBAL_FORCE_RESET_BUTTON, glui_cb);
	global_force_reset_button->set_alignment(GLUI_ALIGN_CENTER);
	glui->add_column_to_panel(panel5, false);
	global_force_direction1_xy_translation =
		glui->add_translation_to_panel(panel5, "direction 1 XY", GLUI_TRANSLATION_XY, &global_force_direction1.x, GLOBAL_FORCE_DIRECTION1_XY_TRANSLATION, glui_cb);
	global_force_direction2_xy_translation =
		glui->add_translation_to_panel(panel5, "direction 2 XY", GLUI_TRANSLATION_XY, &global_force_direction2.x, GLOBAL_FORCE_DIRECTION2_XY_TRANSLATION, glui_cb);
	glui->add_column_to_panel(panel5, false);
	global_force_direction1_z_translation =
		glui->add_translation_to_panel(panel5, "direction 1 Z", GLUI_TRANSLATION_Z, &global_force_direction1.z, GLOBAL_FORCE_DIRECTION1_Z_TRANSLATION, glui_cb);
	global_force_direction2_z_translation =
		glui->add_translation_to_panel(panel5, "direction 2 Z", GLUI_TRANSLATION_Z, &global_force_direction2.z, GLOBAL_FORCE_DIRECTION2_Z_TRANSLATION, glui_cb);

	glui->add_column_to_panel(panel6, false);
	GLUI_Panel *panel3 = glui->add_panel_to_panel(panel6, "Actions/Info");
	screendump_button =
		glui->add_button_to_panel(panel3, "Write images", SCREENDUMP_BUTTON, glui_cb);
	glui->add_button_to_panel(panel3, "Quit", 0, (GLUI_Update_CB)exit);
	glui->add_separator_to_panel(panel3);
	cluster_statictext =
		glui->add_statictext_to_panel(panel3, "");
	cluster_statictext->set_alignment(GLUI_ALIGN_CENTER);
}

// usage message
static void
usage(char *program_name, char *message)
{
	if (message)
		cerr << message << endl;

	cerr << "Usage: " << program_name << endl;
	cerr << "-db\t\tdo double buffering (default)\n";
	cerr << "-sb\t\tdo single buffering\n";
	exit(1);
}

// argument parser
static void
args(int argc, char **argv)
{
	char *program = argv[0];
	
	while (--argc)
	{
		++argv;
		
		if (!strcmp("-h", argv[0]) || !strcmp("-help", argv[0]))
			usage(program, NULL);
		else if (!strcmp("-db", argv[0]) || !strcmp("-double", argv[0]))
			double_buffer = true;
		else if (!strcmp("-sb", argv[0]) || !strcmp("-single", argv[0]))
			double_buffer = false;
		else
			usage(program, "Invalid option!");
	}
}

int
main (int argc, char **argv)
{
	srand48( (unsigned)time( NULL ) );
	glutInit(&argc, argv);
	args(argc, argv);
	
	int type = GLUT_RGB;
#ifdef DEPTH_TEST
	type |= GLUT_DEPTH;
#endif
	type |= GLUT_MULTISAMPLE;
	type |= double_buffer ? GLUT_DOUBLE : GLUT_SINGLE;
	glutInitDisplayMode(type);
	glutInitWindowSize(w, h);
	main_window = glutCreateWindow("COMP 259 Final Project - Max Smolens");

	writer.set_dimensions(w, h);
	global_force_direction1.normalize();
	global_force_direction2.normalize();
	
	glutDisplayFunc(draw);
	glutKeyboardFunc(key_press);
	
	glutCreateMenu(menu);
	glutAddMenuEntry("g: Draw ground", 'g');
	glutAddMenuEntry("k: Draw axes", 'k');
	glutAddMenuEntry("r: Toggle antialiasing", 'r');
	glutAddMenuEntry("p: Use GL_POINTS", 'p');
	glutAddMenuEntry("o: Use textured quad", 'o');
	glutAddMenuEntry("c: Toggle camera motion", 'c');
	glutAddMenuEntry("x: Freeze Particles", 'x');
	glutAddMenuEntry("<esc> or q: exit program", '\033');
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	
#ifdef DEPTH_TEST
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
#endif

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glLineWidth(1.0);
	glPointSize(3.0);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	init_glui();
	step_smoke();
	glutMainLoop();
	return 0;
}
