#ifndef STATE_H
#define STATE_H

#include "stdio.h"

//the state of the current collection of atoms
typedef struct state{
  float N; //how many points
  float *x; //positions
  float E; //potential energy
  int iters;//iterations for basin optimization
}state;

//Local variables, save on allocation time...
namespace State{
  static float statex,statey,statez;
}

void initState(state*,int);
void printState(state*,FILE*);
float dist(float*,float*);
#endif
