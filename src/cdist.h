// cluster distance

#ifndef _CDIST_H
#define _CDIST_H

using namespace std;

class cdist {
public:
	cdist(int _g, float _d) : g(_g), d(_d) { }
	~cdist() {}
	bool operator< (const cdist& c) const {
		return (d < c.d);
	}
	int g;			// particle group id
	float d;		// distance
};

#endif
