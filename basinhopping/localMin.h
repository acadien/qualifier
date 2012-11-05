#ifndef LOCALMIN_H
#define LOCALMIN_H

#include "state.h"

namespace Powell{
  static float **xi;
}

void initPowell(int);
void basinPowell(state* s,float (*func)(float [], void*),void* args);

#endif
