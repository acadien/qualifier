#ifndef LOCALMIN_H
#define LOCALMIN_H

#include "state.h"

void basinPowell(state* s,float ftol, float (*func)(float [], void*),void* args);

#endif
