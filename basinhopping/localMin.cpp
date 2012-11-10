#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "localMin.h"
#include "state.h"
#include "constants.h"
#include "structure.h"

#include "nr.h"
#include "nrutil.h"
#include "powell.c"
#include "linmin.c"
#include "f1dim.c"
#include "mnbrak.c"
#include "brent.c"

using namespace std;

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
  powell(s->x,s->xi,s->N*3,ftol,&(s->iters),&fret,func,args);

  //re-center about the center of mass.
  recenter(s);

  s->E=func(s->x,args);
}
