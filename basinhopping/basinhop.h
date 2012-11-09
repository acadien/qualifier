#ifndef BASINHOP_H
#define BASINHOP_H
#include <queue>
#include "stdlib.h"
#include "state.h"
#include "potential.h"

void MCstep(state*, void*,float,float&,float&,std::queue<int>*,float*,bool);
void resetWindow(std::queue<int>*,float*,int);

#endif
