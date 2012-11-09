#ifndef LOCALMIN_H
#define LOCALMIN_H

#include "state.h"

namespace Powell{
  static float **xi;
  //static const float **xi;
}

void initPowell(int);
void basinPowell(state* s,float ftol, float (*func)(float [], void*),void* args);

#endif
