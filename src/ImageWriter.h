#ifndef _IMAGE_WRITER_H
#define _IMAGE_WRITER_H

#include <iostream>
#include <GL/gl.h>
#include <magick/api.h>

using namespace std;

class ImageWriter
{
public:
	ImageWriter(int w=512, int h=512);
	void output(int frameNum);
	void set_dimensions(int w, int h);
private:
	int w, h;

};

#endif
