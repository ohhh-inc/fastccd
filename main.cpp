// 
// Released by Suejung Huh and Young Ju Lee 2015
//

#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <ctime>
#include <cstdlib>
#include <iostream>
#include "vec4d.h"
#include "ccd_vf.h"


int main()
{
  //ccdTimeResolution = double(1e-7);
  //ccdTimeResolution = double(1e-8);
  ccdTimeResolution = double(1e-9);
  //ccdTimeResolution = double(1e-10);
  //ccdTimeResolution = double(1e-11);

  std::vector <Vec3d> x;
  x.resize(8);
  srand(1);
  for (int j = 0; j < 8; j++ ) {
      x[j].x = ((double)rand())/RAND_MAX;
      x[j].y = ((double)rand())/RAND_MAX;
      x[j].z = ((double)rand())/RAND_MAX;
  }

#ifdef _DEBUG
  special = true; // debug purpose
#endif
   
  //fastccd
  unsigned int numRoots;
  double t[3];
  if (fastccd_Intersect_VF(numRoots,t,x[0], x[1], x[2], x[3],x[4], x[5], x[6], x[7])) {
      std::cout << "collision!\n";
  } else {
      std::cout << "no collision!\n";
  }
}

