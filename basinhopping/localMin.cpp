#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "localMin.h"
#include "state.h"
#include "constants.h"
#include "structure.h"
#include "basinhop.h"
#include "random.h"

#include "nr.h"
#include "nrutil.h"
#include "powell.c"
#include "linmin.c"
#include "f1dim.c"
#include "mnbrak.c"
#include "brent.c"

using namespace std;
using namespace Potential;

//Total system LJ energy
float LJpot(float* cs, void* args){

  float UR=0.0;

  int N=((ARGST*)args)->N;

  for(int i=0;i<N;i++){
    xi=cs[3*i+1];
    yi=cs[3*i+2];
    zi=cs[3*i+3]; 
    for(int j=i+1;j<N;j++){
      xj=cs[3*j+1];
      yj=cs[3*j+2];
      zj=cs[3*j+3];

      dx=xi-xj;
      dy=yi-yj;
      dz=zi-zj;
      //printf("%d %f %f %f\n",j,dx,dy,dz);
      ir=1.0/sqrt(dx*dx+dy*dy+dz*dz);
      ir3=ir*ir*ir;

      UR+=4.0*(ir3*ir3*ir3*ir3 - ir3*ir3);
    }
  }
  return UR;
}

//Single atom's total LJ energy
float LJpotAtom(float* cs, void* args){

  float UR=0.0;

  int N=((ARGST*)args)->N;
  int d=((ARGST*)args)->d;

  xi=cs[3*d+1];
  yi=cs[3*d+2];
  zi=cs[3*d+3]; 
  for(int j=0;j<N;j++){
    xj=cs[3*j+1];
    yj=cs[3*j+2];
    zj=cs[3*j+3];

    dx=xi-xj;
    dy=yi-yj;
    dz=zi-zj;
    ir=1.0/sqrt(dx*dx+dy*dy+dz*dz);
    ir3=ir*ir*ir;

    UR+=4.0*(ir3*ir3*ir3*ir3 - ir3*ir3);
  }
  
  return UR;
}

//If point is outside the sphere, punish the potential
float LJpotPunish(float* cs, void* args){

  float UR=LJpot(cs,args);

  int N=((ARGST*)args)->N;

  for(int i=0;i<N;i++){
    r=origDist(&(cs[3*i+1]));
    
    if( r > boundr )
      UR+=(r-boundr+0.5);
  }

  return UR;
}

//Method A
void basinPowell(state* s,float ftol,float (*func)(float [], void*),void* args){
  float fret;
  powell(s->x,s->xi,s->N*3,ftol,&(s->iters),&fret,func,args);

  //re-center about the center of mass.
  recenter(s);

  s->E=func(s->x,args);
}

void basinJiggle(state* s, float (*func)(float [], void*),void* args){
  //Settings
  int m=s->N*10;
  int k=50;

  //Initialize
  state* sp=((ARGST*)args)->sp;
  state* spp=((ARGST*)args)->spp;
  float alpham=((ARGST*)args)->alpham;
  int N=s->N;
  float oldE;
  s->E=func(s->x,args);

  int tryCount=0,failCount=0;
  float *x=spp->x;

  copyState(s,sp);
  while(failCount < k and tryCount++ < m){
    
    oldE=sp->E;
    
    //Perturb
    for(int i=1;i<=3*N;i++)
      if(mrand()>0.7)
	x[i]=sp->x[i]+(mrand()-0.5)*2*alpham;
    
    spp->E=func(x,args);

    if(oldE > spp->E)
      copyState(spp,sp);
    else
      failCount++;
    printf("*");
  }
  copyState(sp,s);
  if(s->E>1047806272){
    printf("WTF\n");
    exit(0);
  }
  printf("!!! %f\n",alpham);
  recenter(s);
}
