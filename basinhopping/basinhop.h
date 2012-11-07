#ifndef BASINHOP_H
#define BASINHOP_H

#include "state.h"

//Args for powell function call
typedef struct ARGST{
  int N;
}ARGST;

bool MCstep(state*, void*,bool,float);

#endif
