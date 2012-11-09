#ifndef POTENT_H
#define POTENT_H

namespace Potential{
  static float xi,xj,yi,yj,zi,zj;
  static float dx,dy,dz,r;
  static const float sig=1.0;
  static const float feps=4.0*1.0;
}

//Args for LJ/powell function call
typedef struct ARGST{
  int N;
}ARGST;

float LJpot(float* cs, void* args);
float LJpotPunish(float* cs, void* args);

#endif
