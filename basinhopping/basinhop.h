#ifndef BASINHOP_H
#define BASINHOP_H
#include <queue>
#include "stdlib.h"
#include "state.h"
#include "localMin.h"

void MCstep(state*, state*,void*,float,float&,float&,std::queue<int>*,float*,bool);
void resetWindow(std::queue<int>*,float*,int);
void loadIdeal(state* sideal,FILE* ifile);

#endif
