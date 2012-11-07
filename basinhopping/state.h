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

void allocState(state*,int);
void initState(state*,int);
void printState(state*,FILE*);
void printStateBounds(state*, FILE*);
void printStateEnergy(state*, FILE*);
void printStateVolume(state*, FILE*);
void copyState(state*,state*);
#endif
