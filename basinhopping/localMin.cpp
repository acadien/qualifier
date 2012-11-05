#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "localMin.h"

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

void basinPowell(state* s,float (*func)(float [], void*),void* args){
  float fret;
  powell(s->x,xi,s->N*3,1E-4,&(s->iters),&fret,func,args);
  s->E=func(s->x,args);
}
