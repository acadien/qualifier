#include "math.h"
#include "structure.h"
#include "state.h"
#include "constants.h"

using namespace std;
using namespace Struc;

//Distance between 2 points
float dist(float *a, float* b){
  strucx=a[0]-b[0];
  strucy=a[1]-b[1];
  strucz=a[2]-b[2];
  return sqrt(strucx*strucx + strucy*strucy + strucz*strucz);
}

//Distance between point and origin
float origDist(float *a){
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

/*
float salt(state* s){
  int N = s->N;

  for(int i=0;i<N;i++){
    if( origDist(&(s->x[3*i])) > boundr ){
      
    }
  }
  }*/
