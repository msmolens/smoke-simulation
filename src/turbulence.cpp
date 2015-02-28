#include <cstdlib>
#include <cmath>
#include "Noise.h"
#include "turbulence.h"

using namespace std;

float
turbulence(float x, float y, float z)
{
	static Noise n;
	static const int MAX_FREQ = 256;
	int f = 1;
	float turb = 0.0;
	while (f < MAX_FREQ) {
		turb += fabs(1.0/float(f) * (float)n.noise(x, y, z));
		f *= 2;
	}
	return turb;
}
