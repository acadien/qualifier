#ifndef BASINHOP_H
#define BASINHOP_H

#include "state.h"
#include "stdlib.h"
#include <queue>

//Args for powell function call
typedef struct ARGST{
  int N;
}ARGST;

void MCstep(state*, void*,float*,std::queue<int>*,float*);

#endif
