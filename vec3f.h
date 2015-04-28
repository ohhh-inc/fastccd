// 
// Released by Suejung Huh and Young Ju Lee 2015
//

#pragma once

#include <math.h>


#ifdef _DEBUG
#include <stdio.h>
//#include <glh/glh_linear.h>
#endif

#define     GLH_ZERO            float(0.0)
#define     GLH_EPSILON         float(10e-6)
#define		GLH_EPSILON_2		float(10e-12)
#define     equivalent(a,b)     (((a < b + GLH_EPSILON) && (a > b - GLH_EPSILON)) ? true : false)

class Vec3f {
public :
	union {
		struct {
		float x, y, z;
		};
		struct {
		float v[3];
		};
	};
public:

	__inline Vec3f ()
	{x=0; y=0; z=0;}

	__inline Vec3f(const Vec3f &_v)
	{
		x = _v.x;
		y = _v.y;
		z = _v.z;
	}

	__inline Vec3f(const float *_v)
	{
		x = _v[0];
		y = _v[1];
		z = _v[2];
	}

	__inline Vec3f(float _x, float _y, float _z)
	{
		this->x = _x;
		this->y = _y;
		this->z = _z;
	}

	__inline float operator [] ( int i ) const {return v[i];}

	__inline Vec3f &operator += (const Vec3f &_v) {
		x += _v.x;
		y += _v.y;
		z += _v.z;
		return *this;
	}

	__inline Vec3f &operator -= (const Vec3f &_v) {
		x -= _v.x;
		y -= _v.y;
		z -= _v.z;
		return *this;
	}

	__inline void negate() {
		x = -x;
		y = -y;
		z = -z;
	}

	__inline Vec3f operator - () const {
		return Vec3f(-x, -y, -z);
	}

	__inline Vec3f operator+ (const Vec3f &_v) const
	{
		return Vec3f(x+_v.x, y+_v.y, z+_v.z);
	}

	__inline Vec3f operator- (const Vec3f &_v) const
	{
		return Vec3f(x-_v.x, y-_v.y, z-_v.z);
	}

	__inline Vec3f operator *(float t) const
	{
		return Vec3f(x*t, y*t, z*t);
	}

     // cross product
     __inline const Vec3f cross(const Vec3f &vec) const
     {
          return Vec3f(y*vec.z - z*vec.y, z*vec.x - x*vec.z, x*vec.y - y*vec.x);
     }

	 __inline float dot(const Vec3f &vec) const {
		 return x*vec.x+y*vec.y+z*vec.z;
	 }

	 __inline void normalize() 
	 { 
		 float sum = x*x+y*y+z*z;
		 if (sum > GLH_EPSILON_2) {
			 float base = float(1.0/sqrt(sum));
			 x *= base;
			 y *= base;
			 z *= base;
		 }
	 }

	 __inline float length() {
		 return float(sqrt(x*x + y*y + z*z));
	 }

	__inline Vec3f & set_value( const float &vx, const float &vy, const float &vz)
	{ x = vx; y = vy; z = vz; return *this; }

	__inline bool equal_abs(const Vec3f &other) {
		return x == other.x && y == other.y && z == other.z;
	}

	__inline float square_norm() const {
		return x*x+y*y+z*z;
	}

#ifdef _DEBUG
	__inline void output() {
		printf("x=%f, y=%f, z=%f\n", x, y, z);
	}
#endif
};
