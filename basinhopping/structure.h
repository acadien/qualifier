#ifndef STRUCTURE_H
#define STRUCTURE_H

//Local variables, save on allocation time...
#include "state.h"

namespace Struc{
  static float strucx,strucy,strucz;
}

float dist(float *a, float* b);
float origDist(float *a);
float* com(state* s);
float msd(state* a, state* b);
void recenter(state* s);
float salt(state* s);

#endif
