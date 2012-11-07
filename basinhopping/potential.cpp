#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "basinhop.h"
#include "structure.h"
#include "constants.h"
#include "potential.h"

using namespace std;
using namespace Potential;

float LJpot(float* cs, void* args){

  float UR=0.0;

  int N=((ARGST*)args)->N;

  for(int i=0;i<N;i++){
    xi=cs[3*i];
    yi=cs[3*i+1];
    zi=cs[3*i+2]; 
    //printf("%f %f %f\n",xi,yi,zi);
    for(int j=i+1;j<N;j++){
      xj=cs[3*j];
      yj=cs[3*j+1];
      zj=cs[3*j+2];

      dx=xi-xj;
      dy=yi-yj;
      dz=zi-zj;
      r=sqrt(dx*dx+dy*dy+dz*dz);

      UR+=feps*(pow(sig/r,12.0)-pow(sig/r,6.0));
      //break;
    }
    //break;
  }
  
  //printf("atom1:%g %g %g\n",xi,yi,zi);
  //printf("atom2:%g %g %g\n",xj,yj,zj);
  //printf("dist:%g\n",r);
  //printf("energy:%g\n",UR);
 
  return UR;
}

//If point is outside the sphere, punish the potential
float LJpotPunish(float* cs, void* args){

  float UR=LJpot(cs,args);

  int N=((ARGST*)args)->N;
  float r;

  for(int i=0;i<N;i++){
    r=origDist(&(cs[3*i]));
    
    if( r > boundr )
      UR+=(r-boundr+0.5);
  }

  return UR;
}
