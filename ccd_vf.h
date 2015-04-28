// 
// Released by Suejung Huh and Young Ju Lee 2015
//

#pragma once

#include <vector>
#include <set>
#include <algorithm>
#include "vec3f.h"
#include "vec4d.h"
#include <cassert>
#include <cmath>

using namespace std;

//#define NEWTON_RAPHSON
#define HOUSEHOLDER
//#define MACHINE_EPSILON double(1e-15)
//#define MACHINE_EPSILON double(1e-32)
#define MACHINE_EPSILON double(1e-16)

double ccdTimeResolution;
bool special;
int coplanarCase;
int collinearCase;
unsigned int evenCollisions;
int newtonFail;

inline void pr(const Vec3d &x){
  printf("%.25f %.25f %.25f\n",x[0],x[1],x[2]);
}

inline void pr(const string str, const Vec3d &x){
  std::cout << str ;
  printf("%.25f %.25f %.25f\n",x[0],x[1],x[2]);
}

inline Vec3d absAdd(const Vec3d &a, const Vec3d &b) {
  Vec3d res;
  res.x = fabs(a[0]) + fabs( b[0]);
  res.y = fabs(a[1]) + fabs( b[1]);
  res.z = fabs(a[2]) + fabs( b[2]);
  return res;
}

inline double pickMax(const Vec3d &a) {
  if ( a[0] >= a[1] ) {
    if ( a[0] >= a[2] ) {
      return a[0];
    } else 
      return a[2];
  } else {
    if ( a[1] >= a[2] ) {
      return a[1];
    } else 
      return a[2];
  }
}

inline Vec3d absAdd(const Vec3d &a, const Vec3d &b, const Vec3d &c) {
  Vec3d res;
  res.x = fabs(a[0]) + fabs( b[0])+ fabs( c[0]);
  res.y = fabs(a[1]) + fabs( b[1])+ fabs( c[1]);
  res.z = fabs(a[2]) + fabs( b[2])+ fabs( c[2]);
  return res;
}

inline void fabs( Vec3d &a) {
  a.x = fabs(a[0]);
  a.y = fabs(a[1]);
  a.z = fabs(a[2]);
}

#define ISCLOSE_TOLERANCE(a,b,tol) (fabs((a) - (b)) <= tol )
#define ISCLOSE(a,b) (fabs((a) - (b)) <= ccdTimeResolution )

inline bool isClose(const Vec3d &a, const Vec3d &b, const double tol = ccdTimeResolution) {

  if ( !ISCLOSE_TOLERANCE(a[0],b[0],tol)) return false;
  if ( !ISCLOSE_TOLERANCE(a[1],b[1],tol)) return false;
  if ( !ISCLOSE_TOLERANCE(a[2],b[2],tol)) return false;
  return true;
}

inline bool isZero(const Vec3d &a) {
  return ( a[0] == 0 && a[1] == 0 && a[2] == 0 );
}

inline bool isSame(const Vec3d &a, const Vec3d &b) {

  return isZero(a-b);
}

inline bool isNearZero(const Vec3d &a) {
  return ( fabs(a[0]) < ccdTimeResolution && fabs(a[1]) < ccdTimeResolution && fabs(a[2]) < ccdTimeResolution );
}

inline bool isNearZero(const double &a) {
  return ( fabs(a) < ccdTimeResolution );
}

// set bounding box for a point
inline void
setPointBoundingBox(
Vec3d &minAabb, Vec3d &maxAabb,
const Vec3d &q0, const Vec3d &q1)
{
    if ( q0[0] <= q1[0] ) {
      minAabb.x = q0[0];
      maxAabb.x = q1[0];
    } else {
      minAabb.x = q1[0];
      maxAabb.x = q0[0];
    }
    if ( q0[1] <= q1[1] ) {
      minAabb.y = q0[1];
      maxAabb.y = q1[1];
    } else {
      minAabb.y = q1[1];
      maxAabb.y = q0[1];
    }
    if ( q0[2] <= q1[2] ) {
      minAabb.z = q0[2];
      maxAabb.z = q1[2];
    } else {
      minAabb.z = q1[2];
      maxAabb.z = q0[2];
    }

}

// for face and point bounding box test
inline bool  doesFacePointIntersect (const Vec3d &minAabb0, const Vec3d &maxAabb0,
				     const Vec3d &minAabb1, const Vec3d &maxAabb1,
				     const Vec3d &minAabb2, const Vec3d &maxAabb2,
				     const Vec3d &pointMinAabb, const Vec3d &pointMaxAabb)
{
  for ( int i = 0; i < 3 ; i++ ) {
    if ( maxAabb1[i] <= maxAabb0[i] && maxAabb2[i] <= maxAabb0[i] ) {
      if ( maxAabb0[i] + ccdTimeResolution < pointMinAabb[i] ) return false;
    } else if ( maxAabb0[i] <= maxAabb1[i] && maxAabb2[i] <= maxAabb1[i] ) {
      if ( maxAabb1[i] + ccdTimeResolution < pointMinAabb[i] ) return false;
    } else {
      if ( maxAabb2[i] + ccdTimeResolution < pointMinAabb[i] ) return false;
    }

    if ( minAabb0[i] <= minAabb1[i] && minAabb0[i] <= minAabb2[i] ) {
      if ( pointMaxAabb[i] < minAabb0[i] - ccdTimeResolution ) return false;
    } else  if ( minAabb1[i] <= minAabb0[i] && minAabb1[i] <= minAabb2[i] ) {
      if ( pointMaxAabb[i] < minAabb1[i] - ccdTimeResolution ) return false;
    } else {
      if ( pointMaxAabb[i] < minAabb2[i] - ccdTimeResolution ) return false;
    }
  }
  return true;
}


// false : no roots
// true : one or more roots
inline bool solveQuadratic(unsigned int & numRoots, double t[], const double coeffs[])
{
  const double &a = coeffs[2];
  const double &b = coeffs[1];
  const double &c = coeffs[0];

  if ( 0.0 ==a) {

    if ( 0==b ) {
      // root doesn't exist
      numRoots = 0;
      return false;
    }

    t[0]= -c/b;
    numRoots = 1;
    return true;
  }

  double det = b*b - (double)(4.0)* a*c;

  if (det < 0.0 ) {
    // root doesn't exist
    numRoots = 0;
    return false;
  }

  if (det == 0.0 ) {
    // it can have one root
    t[0] = (-b)*0.5/a;
    numRoots = 1;
    return true;
  } 

  // it can have two roots
  det = std::sqrt(det);
  t[0] = (-b - det)*0.5/a;
  t[1] = (-b + det)*0.5/a;
  if ( t[0] > t[1] ) {
    double tmp = t[0];
    t[0] = t[1];
    t[1] = tmp;
  }

  numRoots = 2;
  return true;
}



	
/*
* FASTCCD Areal Condition
*/


static inline bool isInside(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &p, Vec3d &baryc, const double &areaTolerance )
{

	Vec3d n, da, db, dc, dbXdc, dcXda, daXdb;

	Vec3d ba = b-a;
	Vec3d ca = c-a;
	n = ba.cross(ca);
	Vec3d faceArea(fabs(n[0]),fabs(n[1]),fabs(n[2]));

	da = a - p, db = b - p, dc = c - p;
	dbXdc = db.cross(dc);
	dcXda = dc.cross(da);
	daXdb = da.cross(db);
	Vec3d areaA(fabs(dbXdc[0]),fabs(dbXdc[1]),fabs(dbXdc[2]));
	Vec3d areaB(fabs(dcXda[0]),fabs(dcXda[1]),fabs(dcXda[2]));
	Vec3d areaC(fabs(daXdb[0]),fabs(daXdb[1]),fabs(daXdb[2]));

#ifdef _DEBUG
	if ( special ) {
	  printf("area difference = %.30f\n", 
		 faceArea.length() - areaA.length() -areaB.length() -areaC.length());
	}
#endif

	if ( isZero(faceArea) && isZero(areaA) && isZero(areaB) && isZero(areaC) ) {
		// We have to handle degenerated case here too!!!
		Vec3d bc = b-c;
		double baLen = ba.length();
		double caLen = ca.length();
		double bcLen = bc.length();
		Vec3d startPosition, endPosition;
		if ( baLen >= caLen && baLen >= bcLen ) {
			startPosition = a;
			endPosition = b;
		} else if ( caLen >= baLen && caLen >= bcLen ) {
			startPosition = a;
			endPosition = c;
		} else  {
			startPosition = b;
			endPosition = c;
		}
		Vec3d pq = p - startPosition;
		Vec3d lineSegment = endPosition - startPosition;
		if ( pq.dot(lineSegment) < 0 ) {
		  //std::cout << "inside weird 0\n";
		  return false;
		}
		Vec3d crossv = pq.cross(lineSegment);
		if ( crossv.length() != 0 ) {
		  //std::cout << "inside weird 1\n";
		  return false;
		}
		if ( pq.length() > lineSegment.length() ) {
		  //std::cout << "inside weird 2\n";
		  return false;
		}
		return true;
	}

	// regular shape

	double areaAlen = areaA.length();
	double areaBlen = areaB.length();
	double areaClen = areaC.length();
	double faceAreaLength = faceArea.length();

	//if ( !ISCLOSE_TOLERANCE(faceAreaLength,areaAlen+areaBlen+areaClen,ccdTimeResolution)) {
	if ( !ISCLOSE_TOLERANCE(faceAreaLength,areaAlen+areaBlen+areaClen,areaTolerance)) {
#ifdef _DEBUG
	    if ( special ) {
	      std::cout << "inside weird 4 :";
	      printf("%.100f\n",faceAreaLength - ( areaAlen+areaBlen+areaClen));
	    }
#endif
	    return false;
	}

	//Compute barycentric coordinates
	baryc = Vec3d(areaAlen/faceAreaLength, areaBlen/faceAreaLength, areaClen/faceAreaLength);

	return true;
}


// Bisection
inline bool solveCubicBisection(double &l, double &r,const double coeffs[])
{

	// bound the cubic
	double lvalue=coeffs[3]*l*l*l+coeffs[2]*l*l+coeffs[1]*l+coeffs[0];
	double rvalue=coeffs[3]*r*r*r+coeffs[2]*r*r+coeffs[1]*r+coeffs[0];
	if ( lvalue*rvalue > 0 ) {
	  return false;
	}

	// divide and conguer...
	while (fabs(r-l) > ccdTimeResolution) {
		double m=0.5*(r+l);
		double mvalue=coeffs[3]*m*m*m+coeffs[2]*m*m+coeffs[1]*m+coeffs[0];
		if ( lvalue*mvalue < 0 ) {
		  r = m;
		  continue;
		}
	    l = m;
	}
	return true;
}

// Newton-Raphson
inline bool solveCubicNewtonRaphson(double &x,const double coeffs[], const int maxIter)
{

  double deriv,value,change;

  // initial guess x
  int i = 0;
  while ( i < maxIter) {

    deriv = 3.0*coeffs[3]*x*x+2.0*coeffs[2]*x+coeffs[1];
    value = coeffs[3]*x*x*x+coeffs[2]*x*x+coeffs[1]*x+coeffs[0];

    change = -value/deriv;
    x= x + change;

    if ( fabs(change) < ccdTimeResolution ) {
      //if (i > 4 ) std::cout << i << " newton raphson \n";
      return true;
    }
    i++;
  }
  return false;
}

	

// Halley
inline bool solveCubicHouseholder(double &x,const double coeffs[], const int maxIter)
{

 double deriv1,deriv2,value,change;

  // initial guess x
  int i = 0;
  while ( i < maxIter) {

    deriv1 = 3.0*coeffs[3]*x*x+2.0*coeffs[2]*x+coeffs[1];
    deriv2 = 6.0*coeffs[3]*x+2.0*coeffs[2];
    value = coeffs[3]*x*x*x+coeffs[2]*x*x+coeffs[1]*x+coeffs[0];

    change = -(2.0*value*deriv1)/(2.0*deriv1*deriv1-value*deriv2);
    x= x + change;

    if ( fabs(change) < ccdTimeResolution ) {
      //if ( i > 2 ) {
      //      	std::cout << i << " houesholder f(" << x << ")=" << value << ": " << coeffs[3] << "*x*x*x+( " << coeffs[2] << ")*x*x+(" << coeffs[1] << ")*x+" << coeffs[0] << "=0\n";
      //      }
      return true;
    }
    i++;
  }
  return false;
}



inline bool solveCubic(double &l, double &r,const double coeffs[], const int maxIter)
{

  double lvalue=coeffs[3]*l*l*l+coeffs[2]*l*l+coeffs[1]*l+coeffs[0];
  double rvalue=coeffs[3]*r*r*r+coeffs[2]*r*r+coeffs[1]*r+coeffs[0];

  // initial guess
  double x = l + (r-l)*(-lvalue/(rvalue-lvalue));

#ifdef NEWTON_RAPHSON
  if ( solveCubicNewtonRaphson(x,coeffs,maxIter) ) {
    l = x;
    return true;
  };
#endif
#ifdef HOUSEHOLDER
  if ( solveCubicHouseholder(x,coeffs,maxIter) ) {
    l = x;
    return true;
  };
#endif
  return false;
}



// http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline3d/

/*
Calculate the line segment PaPb that is the shortest route between
two lines P1P2 and P3P4. Calculate also the values of mua and mub where
Pa = P1 + mua (P2 - P1)
Pb = P3 + mub (P4 - P3)
Return FALSE if no solution exists.
*/
inline bool LineLineIntersect(
	  const Vec3f &p1, const Vec3f &p2, const Vec3f &p3, const Vec3f &p4,
	  Vec3f &pa, Vec3f &pb, double &mua, double &mub)
{
	Vec3f p13,p43,p21;
	double d1343,d4321,d1321,d4343,d2121;
	double numer,denom;

	p13 = p1 - p3;
	p43 = p4 - p3;
	if (fabs(p43[0])  < GLH_EPSILON && fabs(p43[1])  < GLH_EPSILON && fabs(p43[2])  < GLH_EPSILON)
		return false;

	p21 = p2 - p1;
	if (fabs(p21[0])  < GLH_EPSILON && fabs(p21[1])  < GLH_EPSILON && fabs(p21[2])  < GLH_EPSILON)
		return false;

	d1343 = p13.dot(p43);
	d4321 = p43.dot(p21);
	d1321 = p13.dot(p21);
	d4343 = p43.dot(p43);
	d2121 = p21.dot(p21);

	denom = d2121 * d4343 - d4321 * d4321;
	if (fabs(denom) < GLH_EPSILON)
		return false;
	numer = d1343 * d4321 - d1321 * d4343;

	mua = numer / denom;
	if (mua < 0 || mua > 1)
		return false;

	mub = (d1343 + d4321 * mua) / d4343;
	if (mub < 0 || mub > 1)
		return false;

	pa = p1 + p21*mua;
	pb = p3 + p43*mub;
	return true;
}




// compute a set of intervals which may contain a root
// safeIteration = true if it has a good interval with margin
inline void computeIntervals(const double coeffs[],
			     int &intervalSize,
			   vector < std::pair< double, double > > &intervals)
{

  unsigned int numRoots;
  double t[2];
  double quadCoef[3];
  quadCoef[2] = 3.0*coeffs[3];
  quadCoef[1] = 2.0*coeffs[2];
  quadCoef[0] = coeffs[1];

  solveQuadratic(numRoots,t,quadCoef);

  if ( numRoots == 2 ) {
    if ( ( t[1] <= 0) || (t[0] >= 1 ) || (t[0] <= 0 && t[1] >= 1)) {
      intervalSize = 1;
    } else if ( t[0] <= 0 && 0 <= t[1] && t[1]<= 1 ) {
      intervalSize = 2;
      intervals[0] = make_pair(0.0, t[1]);
      intervals[1] = make_pair(t[1], 1.0);
    } else if (  0 <= t[0] && t[0] <= 1 && 1 <= t[1] ) {
      intervalSize = 2;
      intervals[0] = make_pair(0.0, t[0]);
      intervals[1] = make_pair(t[0], 1.0);
    } else if (  0 <= t[0] && t[0] <= 1 && 0 <= t[1] && t[1] <= 1 ) {
      intervalSize = 3;
      intervals[0] = make_pair(0.0,t[0]);	
      intervals[1] = make_pair(t[0],t[1]);
      intervals[2] = make_pair(t[1],1.0);
    } else {
      std::cout << "fastccd:solveCubic: This can't happen.\n";
      exit(-1);
    }
    return;
  }
   
  //if ( numRoots == 1 ) {
  // numRoots == 0
  // always monotonous no need to partition
  //(  0 <= t[0] && t[0] <= 1 ) 
  //if ( t[0] <= 0  || 1 <=t[0] ) {
  intervalSize = 1;

}

// compute a set of intervals which may contain a root
// safeIteration = true if it has a good interval with margin
inline bool noFlipEarlyOut(const double coeffs[], bool &monotonous)
{

  const double a = 3.0*coeffs[3];
  const double b = 2.0*coeffs[2];
  const double c = coeffs[1];

  const double lvalue = coeffs[0];
  const double rvalue = coeffs[3] + coeffs[2] + coeffs[1] + coeffs[0];
  
  // monotonous
  // and no sign flipping
  // this means there is no zero in [0,1]
  monotonous = ( b*b - 4.0*a*c <= 0 );
  
  if ( monotonous ) {

    if ( lvalue*rvalue > 0 ) return true;

  }

  // there could be zeros
  return false;
}


// root : number of roots
// t : found root time
inline bool solveCubicOneZero(unsigned int & numRoots, double t[], double coeffs[])
{

	numRoots = 0;
	
	double l = 0;
	double r = 1.0;
	double & lvalue = coeffs[0];
	double rvalue = coeffs[3]+coeffs[2]+coeffs[1]+coeffs[0];

	if ( lvalue*rvalue >0 ) return false;

	if (solveCubic(l, r, coeffs,100)) {
	    t[numRoots]= l;
	    numRoots++;
	} else {
	    // ok, newton-raphson failed...
	    newtonFail++;
	    //  std::cout << "Newton-Raphson Failed!!!! Hmmmm!!\n";
	    if (solveCubicBisection(l, r, coeffs)) {
	      t[numRoots]= (l+r)*0.5f;
	      numRoots++;
	    }
	}
	return ( numRoots > 0 );
}



/*
* Computes the coefficients of the cubic equation from the geometry data.
*/
/*
inline void _equateCubic_VF(  const Vec3d &w, const Vec3d &g,
			      const Vec3d &vXu, const Vec3d &vXb_aXu, const Vec3d &aXb,
			      double &a, double &b, double &c, double &d)
{
  a = vXu.dot(w);
  b = vXu.dot(g) + w.dot(vXb_aXu);
  c = aXb.dot(w) + g.dot(vXb_aXu);
  d = aXb.dot(g);
}
*/

inline void _equateCubic_VF(  const Vec3d &A, const Vec3d &B, const Vec3d &C,
			      const Vec3d &U, const Vec3d &V, const Vec3d &W,
			      const Vec3d &bXc, const Vec3d &cXa, const Vec3d &aXb,
			      const Vec3d &vXw, const Vec3d &wXu, const Vec3d &uXv,
			      double &a, double &b, double &c, double &d)
{
  a = vXw.dot(U);
  b = A.dot(vXw) + B.dot(wXu) + C.dot(uXv);
  c = U.dot(bXc) + V.dot(cXa) + W.dot(aXb);
  d = aXb.dot(C);
}


inline void _equateQuad(const Vec3d &alpha, const Vec3d &beta,
			const Vec3d &va, const Vec3d &vb,
			double &a, double &b, double &c)
{
  //Vec3d va, vb;
  //Vec3d alpha, beta;
  Vec3d vaXvb, B, alphaXbeta;

  //va = ad - pd, vb = bd - pd;
  //alpha = a0 - p0, beta = b0 - p0;

  vaXvb = va.cross(vb);
  B = (alpha.cross(vb)+va.cross(beta));
  alphaXbeta = alpha.cross(beta);

#ifdef _DEBUG
  if ( special ) {
    pr(vaXvb);
    pr(B);
    pr(alphaXbeta);
  }
#endif


  if (!(0 ==vaXvb[0] && 0==B[0] && 0==alphaXbeta[0])) {
    a = vaXvb[0], b = -B[0], c = alphaXbeta[0];
  } else if ( !(0==vaXvb[1] && 0==B[1] && 0==alphaXbeta[1])) {
    a = vaXvb[1], b = -B[1], c = alphaXbeta[1];
  } else {
    a = vaXvb[2], b = -B[2], c = alphaXbeta[2];
  } 
}




// compute the time of contact
//
bool
find_contact_time_1d(double &t,
		     const Vec3d &a0,
		     const Vec3d &q0,
		     const Vec3d &qd, const Vec3d &ad)
{
  // Find t where a0 + ad *t = q0 + qd* t
  // t = (q0 - a0)/(ad-qd)

  Vec3d init = q0 - a0;
  Vec3d veld = ad - qd;
  
  if ( isZero(veld) ) return false;

  if ( veld[0] != 0 ) {
    t = init[0]/veld[0];
  } else if ( veld[1] != 0 ) {
    t = init[1]/veld[1];
  } else {
    t = init[2]/veld[2];
  }

#ifdef _DEBUG
  if (special ) {
    std::cout << "1d time " << t << "\n";
  }
#endif

  return ( 0 <= t && t <= 1.0 );
}


// compute the time of contact
//
inline bool
find_contact_2D(unsigned int &numRoots,
		double *t,
		const Vec3d &a0, const Vec3d &b0, const Vec3d &q0, 
		const Vec3d &qd,  const Vec3d &ad,  const Vec3d &bd,
		const Vec3d &oaq, const Vec3d &obq, 
		const Vec3d &daq, const Vec3d &dbq)
{
        //Vec3d qd, ad, bd;
	/* diff. vectors for linear interpolation */
	//qd = q1 - q0, ad = a1 - a0, bd = b1 - b0;

#ifdef _DEBUG
	if ( special ) { std::cout << "\nfind_contact_2D\n"; }
#endif

	/*
	* Compute scalar coefficients by evaluating cross-products.
	*/
	// Form a quadratic equation
	double coeffs[3];
	double a, b, c; /* quadratic polynomial coefficients */
	double tmp;
	double apex, lvalue,avalue,rvalue;

	_equateQuad(oaq, obq, daq, dbq, a, b, c);

	
#ifdef _DEBUG
	if ( special ) { 
	  pr(a0);	  pr(b0);	  pr(q0);
	  std::cout << "coeff=" << a << ", " << b << ", " << c << "\n"; 
	}
#endif

	if (a==0.0 && b==0.0 && c==0.0) {

	  collinearCase++;
#ifdef _DEBUG
	  if ( special ) { std::cout << "1 degenerate\n"; }
#endif
	  
	  // initial inclusion
	  Vec3d ab = a0-b0;
	  fabs(ab);
	  Vec3d lineSeg = absAdd(oaq,obq);
	  if ( isClose(ab,lineSeg) ) {
	    t[0] = 0.0;
	    numRoots = 1;
	    return true;
	  }

	  // dimensional reduction
	  if ( find_contact_time_1d(t[0],a0,q0,qd,ad) ) {
	    if ( find_contact_time_1d(t[1],b0,q0,qd,bd) ) {
	      // two roots case
	      // order
	      if ( t[1] < t[0] ) {
		tmp = t[1];
		t[1] = t[0];
		t[0] = tmp;
	      }
	  
	      // identical case
	      if ( t[0] == t[1] ) {
		numRoots = 1;
		return true;
	      }
	      numRoots = 2;
	      return true;
	    } else {
	      // one root case
	      numRoots = 1;
	      return true;
	    }

	  } else if ( find_contact_time_1d(t[0],b0,q0,qd,bd) ) {

	      // one root case
	      numRoots = 1;
	      return true;
	  }
	  
	  // not intersecting
	  numRoots = 0;
	  return false;
	}

	/*
	 * Well shaped case
	 * Solve the quadratic equation
	 */	
#ifdef _DEBUG
	if ( special ) { std::cout << "Well Shaped 2D ~~~~~\n"; }
#endif
	// check the sign flip
	apex = -0.5*(b/a);
	if ( apex < 0.0 || 1.0 < apex ) {
	  // monotonously increasing or decreasing
	  // [0,1]
	  //double lvalue = c;
	  //double rvalue = a+b+c;
	  if ( c*(a+b+c) > 0 ) {
#ifdef _DEBUG
	    if ( special ) { std::cout << "no sign flip\n"; }
#endif
	    // no sign flip 
	    numRoots = 0;
	    return false;
	  }
	} else {
	  lvalue = c;
	  avalue = a*apex*apex + b*apex + c;
	  rvalue = a+b+c;
	  if ( lvalue*avalue > 0 && avalue*rvalue > 0 ) {
	    // no sign flipping in [0,a] and [a,1]
	    numRoots = 0;
	    return false;
	  }
	}
	
	coeffs[2] = a, coeffs[1] = b, coeffs[0] = c;
	solveQuadratic(numRoots, t, coeffs);

	int validRoots = 0;
	for ( unsigned int i = 0; i< numRoots; i++ ) {
	  if ( 0 <= t[i] && t[i] <= 1.0 ) {
	    
	    t[validRoots] = t[i];
	    validRoots++;
	  }
	}
	numRoots = validRoots;
	return ( numRoots > 0 );
}

inline bool
inPlaneIntersect(const Vec3d &edgeQ,
		 const Vec3d &a0q0, const Vec3d &b0q0, const Vec3d &c0q0, 
		 const Vec3d &a1q0, const Vec3d &b1q0, const Vec3d &c1q0)
{

  Vec3d plane = a0q0.cross(b0q0);
  Vec3d normal = edgeQ.cross(plane);
  double dotA0 = normal.dot(a0q0);

//  if ( isNearZero(plane) ) {
  if ( isZero(plane) ) {
    return true; // can't tell from this info
  }
//  if ( isNearZero(normal) ) {
  if ( isZero(normal) ) {
    return true; // can't tell from this info
  }

  if ( dotA0 > 0 ) {
    if ( normal.dot(a1q0) > 0 && normal.dot(b0q0) > 0 && normal.dot(b1q0) > 0 && normal.dot(c0q0) > 0 && normal.dot(c1q0) > 0 ) return false;
  } else  if ( dotA0 < 0 ) {
    if ( normal.dot(a1q0) < 0 && normal.dot(b0q0) < 0 && normal.dot(b1q0) < 0 && normal.dot(c0q0) < 0 && normal.dot(c1q0) < 0 ) return false;
  }
      
  return true;
}

// isInside between two points
inline bool
isInside(const Vec3d &a,  const Vec3d &b,  const Vec3d &q) {

  for ( int i = 0; i < 3; i++ ) {
    if ( a[i] < b[i] ) {
      if ( q[i] < a[i] || b[i] < q[i] ) return false;
    } else {
      if ( q[i] < b[i] || a[i] < q[i] ) return false;
    }
  }
  Vec3d aq = a - q;
  Vec3d qb = q - b;

  if ( aq[0] != 0 && qb[0] != 0 ) {
    if ( aq[1]/aq[0] != qb[1]/qb[0] ) return false;
    if ( aq[2]/aq[0] != qb[2]/qb[0] ) return false;
    return true;
  } else if ( aq[1] != 0 && qb[1] != 0 ) {
    if ( aq[0]/aq[1] != qb[0]/qb[1] ) return false;
    if ( aq[2]/aq[1] != qb[2]/qb[1] ) return false;
    return true;
  } else if ( aq[2] != 0 && qb[2] != 0 ) {
    if ( aq[0]/aq[2] != qb[0]/qb[2] ) return false;
    if ( aq[1]/aq[2] != qb[1]/qb[2] ) return false;
    return true;
  } 

  return true;
}


inline bool
find_contact_2D(unsigned int& numRoots,
		double *t,
		const Vec3d &a0,  const Vec3d &b0,  const Vec3d &c0,  const Vec3d &q0,
		const Vec3d &a1,  const Vec3d &b1,  const Vec3d &c1,  
		const Vec3d &ad,  const Vec3d &bd,  const Vec3d &cd,  const Vec3d &qd,
		const Vec3d &oaq, const Vec3d &obq, const Vec3d &ocq, 
		const Vec3d &daq, const Vec3d &dbq, const Vec3d &dcq)
{

  numRoots = 0;

  Vec3d a1q0 = a1-q0;
  Vec3d b1q0 = b1-q0;
  Vec3d c1q0 = c1-q0;

  if ( !inPlaneIntersect(qd,oaq,obq,ocq,a1q0,b1q0,c1q0) ) {
    return false;
  }

  std::vector <double> roots;
  roots.resize(6,0.0);
  int rootCount = 0;
  Vec3d a,b,q;

  // dimensional reduction
  if ( find_contact_2D(numRoots, t, a0, b0, q0, qd,ad,bd,oaq,obq,daq,dbq)) {
#ifdef _DEBUG
      if ( special ) std::cout << "case0 collision C " << numRoots << " roots\n";
#endif
    for (unsigned int i = 0; i < numRoots; i++) {
      if ( isInside(a0+ad*t[i],b0+bd*t[i],q0+qd*t[i])) {
#ifdef _DEBUG
	  if ( special ) {
	    std::cout << rootCount << "th root " << t[i] <<"\n";
	    std::cout << "a="; pr(a0+ad*t[i]);
	    std::cout << "b="; pr(b0+bd*t[i]);
	    std::cout << "p="; pr(q0+qd*t[i]);
	  }
#endif
	roots[rootCount] = t[i];
	rootCount++;
      }
    }
  }

  if (find_contact_2D(numRoots, t, b0, c0, q0, qd,bd,cd,obq,ocq,dbq,dcq)) {
#ifdef _DEBUG
      if ( special ) std::cout << "case1 collision C " << numRoots << " roots\n";
#endif
    for (unsigned int i = 0; i < numRoots; i++) {
	if ( isInside(b0+bd*t[i],c0+cd*t[i],q0+qd*t[i])) {
#ifdef _DEBUG
	  if ( special ) {
	    std::cout << rootCount << "th root " << t[i] <<"\n";
	    std::cout << "b="; pr(b0+bd*t[i]);
	    std::cout << "c="; pr(c0+cd*t[i]);
	    std::cout << "p="; pr(q0+qd*t[i]);
	  }
#endif
	  roots[rootCount] = t[i];
	  rootCount++;
	}
      }
  };

  if (find_contact_2D(numRoots, t, c0, a0, q0, qd,cd,ad,ocq,oaq,dcq,daq)) {
#ifdef _DEBUG
      if ( special ) std::cout << "case2 collision C " << numRoots << " roots\n";
#endif

      for (unsigned int i = 0; i < numRoots; i++) {
	if ( isInside(a0+ad*t[i],c0+cd*t[i],q0+qd*t[i])) {
#ifdef _DEBUG
	  if ( special ) {
	    std::cout << rootCount << "th root " << t[i] <<"\n";
	    std::cout << "a="; pr(a0+ad*t[i]);
	    std::cout << "c="; pr(c0+cd*t[i]);
	    std::cout << "p="; pr(q0+qd*t[i]);
	  }
#endif
	  roots[rootCount] = t[i];
	  rootCount++;
	}
      }
  }

  numRoots = rootCount;
  if (numRoots == 0) {
#ifdef _DEBUG
      if ( special ) std::cout << "NO COLLISION\n";
#endif
    return false;
  }

  // Finding too many roots... This is not an error
  //if ( roots.size() > 3 ) {
  //numRoots = roots.size();
  //std::cout << "find_contact_2d : Find too many roots!\n";
  //}

  if ( numRoots == 1 ) {
	t[0] = roots[0];
  } else {
    
    roots.resize(numRoots);
    std::sort(roots.begin(),roots.end());
    rootCount = 1;
    t[0] = roots[0];
    for (unsigned int i = 1; i < numRoots; i++) {
      if ( t[rootCount-1] >= roots[i] ) continue;
      t[rootCount] = roots[i];
      rootCount++;
    }
    numRoots = rootCount;
  }
#ifdef _DEBUG
  if ( special ) {
    std::cout << "numRoots = " << numRoots << "\n";
  }
#endif

  return (numRoots > 0);
}

bool inclusionTest(unsigned int & numContact,
		   double contactTime[],
		   const unsigned int & numRoots,
		   const double t[],
		   const Vec3d &q0, const Vec3d &a0, const Vec3d &b0, const Vec3d &c0,
		   const Vec3d &qd, const Vec3d &ad, const Vec3d &bd, const Vec3d &cd) {

  double y0x0 = (b0-a0).length();
  double z0x0 = (c0-a0).length();
  double vyvx = (bd-ad).length();
  double vzvx = (cd-ad).length();

  double aq0 = (a0-q0).length();
  double bq0 = (b0-q0).length();
  double cq0 = (c0-q0).length();
  double aqv = (ad-qd).length();
  double bqv = (bd-qd).length();
  double cqv = (cd-qd).length();

  double areaTolerance = 0.5*ccdTimeResolution
	  *(y0x0*vzvx+z0x0*vyvx+2.0*vyvx*vzvx
	    +aq0*bqv+bq0*aqv+2.0*aqv*bqv
	    +bq0*cqv+cq0*bqv+2.0*bqv*cqv
	    +cq0*aqv+aq0*cqv+2.0*cqv*aqv)
           +ccdTimeResolution*ccdTimeResolution*(vyvx*vzvx + aqv*bqv + bqv*cqv + cqv*aqv);

  if ( areaTolerance < ccdTimeResolution ) areaTolerance = ccdTimeResolution;

  for ( unsigned int i = 0; i < numRoots; i++ ) 
    {
	if ( t[i] < 0.0 || 1.0 < t[i] ) continue;  

#ifdef _DEBUG
	if ( special ) {
	  std::cout << i << "th " << t[i] << "\n";
	}
#endif
	
	Vec3d ta = a0 + ad*t[i];
	Vec3d tb = b0 + bd*t[i];
	Vec3d tc = c0 + cd*t[i];
	Vec3d q = q0 + qd*t[i];
	Vec3d baryc;

#ifdef _DEBUG
	if ( special ) {
	  printf("areaTolerance = %.30f\n", areaTolerance);
	}
#endif

	if ( isInside(ta,tb,tc,q,baryc,areaTolerance) ) {

#ifdef _DEBUG
	if ( special ) {
	  std::cout << "COLLISION\n";
	}
#endif
	  contactTime[numContact] = t[i];
	  numContact++;
#ifdef _DEBUG
	} else {
	  if ( special ) {
	    std::cout << "NO COLLISION\n";
	  }
#endif
	}
    }
	  
 if ( numContact == 2  ) {
    evenCollisions++;
    return true;
  }
  return (numContact > 0);
}

bool degenerateTest(double contactTime[],
		    unsigned int & numContact,
		    unsigned int & numRoots,
		    double t[],
		    const Vec3d &q0, const Vec3d &a0, const Vec3d &b0, const Vec3d &c0,
		    const Vec3d &a1, const Vec3d &b1, const Vec3d &c1,
		    const Vec3d &qd, const Vec3d &ad, const Vec3d &bd, const Vec3d &cd,
		    const Vec3d &oaq, const Vec3d &obq, const Vec3d &ocq, 
		    const Vec3d &daq, const Vec3d &dbq, const Vec3d &dcq,
		    const Vec3d &obqXocq)
		    
{

   // coplanar
    // degenerate case
    // There could be more than one contact, 
    // However in degenerate case, we only find the earliest contact
    Vec3d oaqXobq = oaq.cross(obq);
    Vec3d ocqXoaq = ocq.cross(oaq);
    Vec3d area = absAdd(oaqXobq,obqXocq,ocqXoaq);
    Vec3d areaTri = (b0-a0).cross(c0-a0);
    fabs(areaTri);

#ifdef _DEBUG
    if (special) {
      pr("areaTri", areaTri);
      pr("area", area);
      pr("oaq",oaq);
      pr("obq",obq);
      pr("ocq",ocq);
      pr("oaqXobq",oaqXobq);
      pr("obqXocq",obqXocq);
      pr("ocqXoaq",ocqXoaq);
      double area00 = sqrt(areaTri[0]*areaTri[0] + areaTri[1]*areaTri[1] + areaTri[2]*areaTri[2]);
      double area11 = sqrt(area[0]*area[0] + area[1]*area[1] + area[2]*area[2]);
      std::cout << area00 << "," << area11 << "\n";
    }
#endif
  
    if ( isClose(areaTri,area) && !isZero(area)) {

      contactTime[0] = 0;
      numContact = 1;
#ifdef _DEBUG
	  if ( special ) { std::cout << "============================== FOUND COLLISION at 0\n"; }
#endif

      return true;
    } else {
#ifdef _DEBUG
	  if ( special ) { 
	    std::cout << "============================== NO COLLISION at 0\n"; 
	  }
#endif

    }

    numContact = 0;
    if ( find_contact_2D(numRoots,t, a0,b0,c0,q0,a1,b1,c1,ad,bd,cd,qd,
			 oaq,obq,ocq,daq,dbq,dcq) ) {

      for ( unsigned int i = 0; i < numRoots; i++ ) 
	{

	  if ( t[i] < 0.0 || 1.0 < t[i] ) continue;  

#ifdef _DEBUG
	  if ( special ) { std::cout << "============================== FOUND COLLISION at " << t[i] << "\n"; }
#endif
	  contactTime[numContact] = t[i];
	  numContact++;

	}
	  
#ifdef _DEBUG
	  if ( special ) { 
	    std::cout << "============================== FOUND " << numContact << " COLLISION \n"; 
	  }
#endif
      // in case of boundary collision, count all as collision cases
      return (numContact > 0);
    }
    return false;
}


// detect contact in 3D and the time of contact
//

bool
find_contact_3D(unsigned int & numContact,
		double contactTime[],
		const Vec3d &q0, const Vec3d &a0, const Vec3d &b0, const Vec3d &c0,
		const Vec3d &q1, const Vec3d &a1, const Vec3d &b1, const Vec3d &c1)
{

  Vec3d qd, ad, bd, cd;
  double t[20];
  unsigned int numRoots;
  bool checkRounding_quadratic = false;
  bool checkRounding_linear = false;
  bool checkRounding_degenerate = false;

  if ( isSame(q0,a0) || isSame(q0,b0) || isSame(q0,c0) ) {
    numContact = 1;
    contactTime[0] = 0;
    return true;
  } else if (  isSame(q1,a1) || isSame(q1,b1) || isSame(q1,c1) ) {
    numContact = 1;
    contactTime[0] = 1.0;
    return true;
  }
    

  /* diff. vectors for linear interpolation */
  qd = q1 - q0, ad = a1 - a0, bd = b1 - b0, cd = c1 - c0;
  numContact = 0;

#ifdef _DEBUG
  if ( special ) {
    std:: cout << "....................................................................velocity\n";
    printf("%.30f %.30f %.30f\n",qd[0],qd[1],qd[2]);
    printf("%.30f %.30f %.30f\n",ad[0],ad[1],ad[2]);
    printf("%.30f %.30f %.30f\n",bd[0],bd[1],bd[2]);
    printf("%.30f %.30f %.30f\n",cd[0],cd[1],cd[2]);
    std::cout << endl ;

    Vec3d ba,ca,p0a,p1a, normal,aq,bq,cq,tri0,tri1,tri2;
    double p0dot, p1dot;
    ba = b0-a0; ca = c0-a0; p0a = q0 - a0; normal = ba.cross(ca); 
    std::cout << "Triangle area at time 0 =" << normal.length() << "\n";
    aq = a0 - q0; bq = b0 - q0; cq = c0 - q0;
    tri0= aq.cross(bq); tri1 = bq.cross(cq); tri2 = cq.cross(aq);
    std::cout << "3 Triangle area at time 0 =" << tri0.length() + tri1.length() + tri2.length() << "\n";
    double areaDiff = tri0.length() + tri1.length() + tri2.length() - normal.length();
    printf(" areaDiff = %.30f\n",areaDiff);
    normal.normalize(); 
    p0dot = p0a.dot(normal);
    printf(" p0dot = %.30f\n",p0dot);
    printf("%.30f %.30f %.30f\n",normal[0],normal[1],normal[2]);

    ba = b1-a1; ca = c1-a1; p1a = q1 - a1; normal = ba.cross(ca); 
    std::cout << "Triangle area at time 1 =" << normal.length();
    aq = a1 - q1; bq = b1 - q1; cq = c1 - q1;
    tri0= aq.cross(bq); tri1 = bq.cross(cq); tri2 = cq.cross(aq);
    std::cout << "3 Triangle area at time 1 =" << tri0.length() + tri1.length() + tri2.length() << "\n";
    areaDiff = tri0.length() + tri1.length() + tri2.length() - normal.length();
    printf(" areaDiff = %.30f\n",areaDiff);
    normal.normalize(); 
    p1dot = p1a.dot(normal);
    printf(" p1dot = %.30f\n",p1dot);
    printf("%.30f %.30f %.30f\n",normal[0],normal[1],normal[2]);
    std::cout << endl << endl;
    
  }
#endif

  /*
   * Compute scalar coefficients by evaluating dot and cross-products.
   */
  // Form a cubic equation

  Vec3d oaq, obq, ocq;
  Vec3d daq, dbq, dcq;
  //Vec3d dbqXdcq, dbqXocq, obqXdcq, obqXocq;
  //Vec3d dbqXocq_obqXdcq;
  Vec3d dbqXdcq, dcqXdaq, daqXdbq;
  Vec3d obqXocq, ocqXoaq, oaqXobq;

  double coeffs[4];
  bool foundRoot = false;
  int intervalSize;
  vector < std::pair< double, double > > intervals;
  intervals.resize(3);
  
  oaq = a0 - q0, obq = b0 - q0, ocq = c0 - q0;
  daq = ad - qd, dbq = bd - qd, dcq = cd - qd;

  dbqXdcq = dbq.cross(dcq);
  dcqXdaq = dcq.cross(daq);
  daqXdbq = daq.cross(dbq);

  obqXocq = obq.cross(ocq);
  ocqXoaq = ocq.cross(oaq);
  oaqXobq = oaq.cross(obq);

  double a, b, c, d; /* cubic polynomial coefficients */

#ifdef _DEBUG
  if ( special ) {
    Vec3d dcqXdaq, daqXdbq;
    Vec3d ocqXoaq, oaqXobq;
    std:: cout << "......................................................... internal values\n";

    dcqXdaq = dcq.cross(daq);
    daqXdbq = daq.cross(dbq);
    ocqXoaq = ocq.cross(oaq);
    oaqXobq = oaq.cross(obq);

    pr("a0-q0           ",oaq);
    pr("b0-q0           ",obq);
    pr("c0-q0           ",ocq);
    pr("V_a-V_q         ",daq);
    pr("V_b-V_q         ",dbq);
    pr("V_c-V_q         ",dcq);
    pr("(b0-q0)x(c0-q0) ",obqXocq);
    pr("(c0-q0)x(a0-q0) ",ocqXoaq);
    pr("(a0-q0)x(b0-q0) ",oaqXobq);
    pr("(Vb-Vq)x(Vc-Vq) ",dbqXdcq);
    pr("(Vc-Vq)x(Va-Vq) ",dcqXdaq);
    pr("(Va-Vq)x(Vc-Vq) ",daqXdbq);
    std::cout << "((b0-q0)x(c0-q0)).dot(Va-Vq)= " << (obqXocq).dot(daq) << "\n";
    std::cout << "((c0-q0)x(a0-q0)).dot(Vb-Vq)= " << (ocqXoaq).dot(dbq) << "\n";
    std::cout << "((a0-q0)x(b0-q0)).dot(Vc-Vq)= " << (oaqXobq).dot(dcq) << "\n";
    std::cout << "((Vb-Vq)x(Vc-Vq)).dot(a0-q0)= " << (dbqXdcq).dot(oaq) << "\n";
    std::cout << "((Vc-Vq)x(Va-Vq)).dot(b0-q0)= " << (dcqXdaq).dot(obq) << "\n";
    std::cout << "((Va-Vq)x(Vb-Vq)).dot(c0-q0)= " << (daqXdbq).dot(ocq) << "\n";


  }
#endif

  _equateCubic_VF(oaq,obq,ocq,daq,dbq,dcq,
		  obqXocq, ocqXoaq,oaqXobq, 
		  dbqXdcq,dcqXdaq,daqXdbq,
		  a, b, c, d);


  coeffs[3] = a, coeffs[2] = b, coeffs[1] = c, coeffs[0] = d;

#ifdef _DEBUG
    if ( special ) {
      std::cout << "A=" << a << " B=" << b << "C=" << c << "D=" << d << "\n";
    }
#endif

  //
  // Is this coplanar case?
  //
  if (a ==0.0 && b==0.0 && c==0.0 && d == 0.0) {

#ifdef _DEBUG
    if ( special ) {
      std::cout << "Coplanar\n"; 
    }
#endif

    coplanarCase++;
 
    return degenerateTest(contactTime,
			  numContact,
			  numRoots,t, q0,a0,b0,c0,a1,b1,c1,qd,ad,bd,cd,
			  oaq,obq,ocq,daq,dbq,dcq, obqXocq);
  } else {

    
    //
    // A well shaped case, i.e. not coplanar
    //

#ifdef _DEBUG
    if ( special ) {
      std::cout << "Not Coplanar\n"; 
    }
#endif

    if ( a == 0.0 ) {
      
      if ( -MACHINE_EPSILON <= b && b <= MACHINE_EPSILON && b != 0 ) {
	if ( c != 0 ) checkRounding_linear = true;
      }

      foundRoot = solveQuadratic(numRoots,t,coeffs);
      
#ifdef _DEBUG
      if ( special ) {
	std::cout << "zero A\n"; 
	if ( foundRoot ) {
	  std::cout << "found " << numRoots << " root \n";
	  printf("t[0] = %.50f\n",t[0]);
	  printf("t[1] = %.50f\n",t[1]);
	} else {
	  std::cout << "found no root\n";
	}
      }
#endif
      
    } else {

      if ( -MACHINE_EPSILON <= a && a <= MACHINE_EPSILON ) {
	if ( -MACHINE_EPSILON <= b && b <= MACHINE_EPSILON ) {
	  if ( -MACHINE_EPSILON <= c && c <= MACHINE_EPSILON ) {
	    if ( -MACHINE_EPSILON <= d && d <= MACHINE_EPSILON ) { 
	      checkRounding_degenerate = true;
	    }
	  } else {
	    if ( c != 0 ) checkRounding_linear = true;
	  }
	} else {
	  checkRounding_quadratic = true;
	}
      }

      // 
      // Well Shpaed Case with A with significant size
      // Solve the cubic equation
      // iteratively solve the cubic (scalar) equation and test for validity of the solution.

#ifdef _DEBUG
      if ( special ) {
	std::cout << "Well Shaped\n"; 
	std::cout << "A=" << a << " B=" << b << "C=" << c << "D=" << d << "\n";
      }
#endif
      double lvalue=coeffs[0];
      double rvalue=coeffs[3]+coeffs[2]+coeffs[1]+coeffs[0];

      bool monotonous = false;
      if ( !(noFlipEarlyOut(coeffs, monotonous)) ) 
      {
      
	if ( monotonous) {
	  foundRoot = solveCubicOneZero(numRoots, t, coeffs);

	} else {
	
	  numRoots = 0;
	  computeIntervals(coeffs,intervalSize, intervals);
	  if (intervalSize == 1) {
	    // at most one root
	    foundRoot = solveCubicOneZero(numRoots, t, coeffs);

	  } else if ( intervalSize == 2 ) {

	    // could have up to three roots
	
	    double l = 0.0;
	    double m = intervals[0].second;
	    double r = 1.0;

	    lvalue=coeffs[0];
	    double mvalue=coeffs[3]*m*m*m+coeffs[2]*m*m+coeffs[1]*m+coeffs[0];
	    rvalue=coeffs[3]+coeffs[2]+coeffs[1]+coeffs[0];
	
	    if ( lvalue * mvalue <= 0 ) {
	      foundRoot = true;
	      if (solveCubicBisection(l, m, coeffs)) {
		t[numRoots]= (l+m)*0.5f;
		numRoots++;
	      }
	    } 
	    if ( mvalue *rvalue <= 0 ) {
	      foundRoot = true;
	      if (solveCubicBisection(m, r, coeffs)) {
		t[numRoots]= (m+r)*0.5f;
		numRoots++;
	      }
	    } 

	  } else {

	    // intervalSize == 3
	    // could have up to three roots

	    double l = 0.0;
	    double m0 = intervals[0].second;
	    double m1 = intervals[1].second;
	    double r = 1.0;

	    lvalue=coeffs[0];
	    double m0value=coeffs[3]*m0*m0*m0+coeffs[2]*m0*m0+coeffs[1]*m0+coeffs[0];
	    double m1value=coeffs[3]*m1*m1*m1+coeffs[2]*m1*m1+coeffs[1]*m1+coeffs[0];
	    rvalue=coeffs[3]+coeffs[2]+coeffs[1]+coeffs[0];
		
	    if ( lvalue * m0value <= 0 ) {
	      foundRoot = true;
	      if (solveCubicBisection(l, m0, coeffs)) {
		t[numRoots]= (l+m0)*0.5f;
		numRoots++;
	      }
	    } 
	    if ( m0value *m1value <= 0 ) {
	      foundRoot = true;
	      if (solveCubicBisection(m0, m1, coeffs)) {
		t[numRoots]= (m0+m1)*0.5f;
		numRoots++;
	      }
	    } 

	    if ( m1value *rvalue <= 0 ) {
	      foundRoot = true;
	      if (solveCubicBisection(m1, r, coeffs)) {
		t[numRoots]= (m1+r)*0.5f;
		numRoots++;
	      }
	    } 
	  }
	}
      }
    }
    if ( foundRoot > 0 ) {
      if ( inclusionTest(numContact,contactTime,
			 numRoots, t,
			 q0,a0,b0,c0,qd,ad,bd,cd) ) return true;
    }
  }

#ifdef HANDLE_ROUNDING

  if ( !checkRounding_degenerate && !checkRounding_linear && !checkRounding_quadratic )  return false;

  // Check if the significant coefficients are too small
  // that they could be non-zeros due to accumulated rounding errors
  // when they are actually zeros.

  if ( checkRounding_linear ) {
#ifdef _DEBUG	  
      if ( special ) {
	std :: cout << "FASTCCD ================================ linear \n";
      }
#endif
      t[0] = -d/c;
      if ( 0<= t[0] && t[0] <= 1 ) {
	numRoots = 1;
	return ( inclusionTest(numContact,contactTime,
			       numRoots, t,
			       q0,a0,b0,c0,qd,ad,bd,cd));
      }
  }

  if ( checkRounding_quadratic ) {
#ifdef _DEBUG	  
    if ( special ) {
      std :: cout << "FASTCCD ============================== quadratic \n";
    }
#endif
    coeffs[3] = 0;
    foundRoot = solveQuadratic(numRoots,t,coeffs);
    return inclusionTest(numContact,contactTime,
			 numRoots, t,
			 q0,a0,b0,c0,qd,ad,bd,cd);
  }

  if ( checkRounding_degenerate ) {
    // test for degenerate case
#ifdef _DEBUG	  
    if ( special ) {
      std :: cout << coeffs[1]  << " FASTCCD ================================ degenerate \n";
    }
#endif
    return degenerateTest(contactTime,
			  numContact,
			  numRoots,t, q0,a0,b0,c0,a1,b1,c1,qd,ad,bd,cd,
			  oaq,obq,ocq,daq,dbq,dcq, obqXocq);

  };

#endif


  return false;
    
}


double
fastccd_Intersect_VF(unsigned int numRoots,
		     double t[],
		     const Vec3d &q0, const Vec3d &a0, const Vec3d &b0, const Vec3d &c0, 
		     const Vec3d &q1, const Vec3d &a1, const Vec3d &b1, const Vec3d &c1)
{

  // Detect a contact for [0,1]
  find_contact_3D(numRoots,t,q0,a0,b0,c0,q1,a1,b1,c1);

#ifdef IGNORE_EVEN
  if ( numRoots == 2 ) return false;
#endif

  return ( numRoots > 0 );
}

