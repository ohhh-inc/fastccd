// 
// Released by Suejung Huh and Young Ju Lee 2015
//

#pragma once

#include <assert.h>
#include <cmath>


class Vec3d {
public:
	union {
		struct {
		double x, y, z;
		};
		struct {
		double v[3];
		};
	};

	__inline Vec3d ()
	{x=0; y=0; z=0;}

	__inline Vec3d(const Vec3d &_v)
	{
		x = _v[0];
		y = _v[1];
		z = _v[2];
	}

	__inline const double & operator [] ( int i ) const
	{
		switch (i) {
		case 0: return x;
		case 1: return y;
		case 2: return z;
		default: assert(0); return x;
		}
	}

	__inline Vec3d(double _x, double _y, double _z)
	{
		this->x = _x;
		this->y = _y;
		this->z = _z;
	}

	__inline Vec3d operator+ (const Vec3d &_v) const
	{
		return Vec3d(x+_v.x, y+_v.y, z+_v.z);
	}

	__inline Vec3d operator- (const Vec3d &_v) const
	{
		return Vec3d(x-_v.x, y-_v.y, z-_v.z);
	}

	__inline Vec3d operator *(double t) const
	{
		return Vec3d(x*t, y*t, z*t);
	}

	__inline bool operator ==(const Vec3d &_v) const
	{
	  return ( x == _v.x && y == _v.y && z == _v.z );
	}

     // cross product
     __inline const Vec3d cross(const Vec3d &vec) const
     {
          return Vec3d(y*vec.z - z*vec.y, z*vec.x - x*vec.z, x*vec.y - y*vec.x);
     }

	 __inline double dot(const Vec3d &vec) const {
		 return x*vec.x+y*vec.y+z*vec.z;
	 }

	 __inline void normalize() 
	 { 
		 double sum = x*x+y*y+z*z;
		 if (sum > 0) {
			 double base = double(1.0/sqrt(sum));
			 x *= base;
			 y *= base;
			 z *= base;
		 }
	 }

	 __inline double length() {
		 return sqrt(x*x + y*y + z*z);
	 }

#ifdef _DEBUG
	__inline void output() const {
		printf("x=%.30f, y=%.30f, z=%.30f\n", (double)x, (double)y, (double)z);
	}
#endif
};
