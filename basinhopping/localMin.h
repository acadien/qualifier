#ifndef LOCALMIN_H
#define LOCALMIN_H

#include "state.h"

namespace Potential{
  static float xi,xj,yi,yj,zi,zj;
  static float dx,dy,dz,r,ir,ir3;
}

//Args for LJ/powell function call
typedef struct ARGST{
  int N;
}ARGST;

float LJpot(float* cs, void* args);
float LJpotPunish(float* cs, void* args);

void basinPowell(state* s,float ftol, float (*func)(float [], void*),void* args);
void basinJiggle(state* s, float ftol, float (*func)(float [], void*),void* args);

#endif
