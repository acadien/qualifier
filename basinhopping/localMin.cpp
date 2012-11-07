#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "localMin.h"
#include "state.h"
#include "constants.h"

#include "nr.h"
#include "nrutil.h"
#include "powell.c"
#include "linmin.c"
#include "f1dim.c"
#include "mnbrak.c"
#include "brent.c"

using namespace std;
using namespace Powell;

void initPowell(int NDIM){
  int i,j;

  xi=matrix(1,NDIM,1,NDIM);  
  for (i=1;i<=NDIM;i++){
    for (j=1;j<=NDIM;j++){
      xi[i][j]=(i==j ? 1.0 : 0.0);
    }
  }
}

/* The original Powell's method, unbounded
void basinPowell(state* s,float ftol,float (*func)(float [], void*),void* args){
  float fret;
  powell(s->x,xi,s->N*3,ftol,&(s->iters),&fret,func,args);
  //  printf("%d powell steps\n",s->iters);
  s->E=func(s->x,args);
}
*/


//Bounded powell's... kind of
void basinPowell(state* s,float ftol,float (*func)(float [], void*),void* args){
  float fret;
  float origin[3]={0.,0.,0.};
  powell(s->x,xi,s->N*3,ftol,&(s->iters),&fret,func,args);

  //First re-center about the center of mass.
  float com[3]={0.,0.,0.};
  for(int i;i<s->N;i++){
    com[0]+=s->x[3*i];
    com[1]+=s->x[3*i+1];
    com[2]+=s->x[3*i+2];
  }
  for(int i;i<3;i++)
    com[i]/=s->N;
  for(int i;i<s->N;i++){
    s->x[3*i]-=com[0];
    s->x[3*i+1]-=com[1];
    s->x[3*i+2]-=com[2];
  }

  s->E=func(s->x,args);
  //How many atoms outside the sphere?
  //int count=0;
  
    
    //if(dist((float*)origin,&(s->x[3*i])) > boundr){
    //  count++;
    //}
  
  //printf("outside atoms %d\n",count);
}
