#ifndef POTENT_H
#define POTENT_H

namespace Potential{
  static float xi,xj,yi,yj,zi,zj;
  static float dx,dy,dz,r;
  static const float sig=1.0;
  static const float feps=4.0*1.0;
}
float LJpot(float* cs, void* args);

#endif
