#include "math.h"
#include "structure.h"
#include "state.h"
#include "constants.h"
#include "potential.h"
#include "random.h"

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

//Calculate the com
float* com(state* s){
  //First re-center about the center of mass.
  float* com=s->com;
  for(int i=0;i<s->N;i++){
    com[0]+=s->x[3*i];
    com[1]+=s->x[3*i+1];
    com[2]+=s->x[3*i+2];
  }
  for(int i=0;i<3;i++)
    com[i]/=s->N;
  return (float*)com;
}

void recenter(state* s){
  float *com=s->com;
  for(int i=0;i<s->N;i++){
    s->x[3*i]-=com[0];
    s->x[3*i+1]-=com[1];
    s->x[3*i+2]-=com[2];
  }
  for(int i=0;i<3;i++)
    s->com[i]=0.;
}

float msd(state* a, state* b){
  float d=0,m=0;
  float* xa=a->x;
  float* xb=b->x;
  int N=a->N;
  for(int i=0;i<N;i++){
    d=dist(&xa[3*i],&xb[3*i]);
    m+=d*d;
  }
  return sqrt(m/(float)N);
}

float salt(state* s){
  int N = s->N,cnt,outcount=0;
  float *x=s->x;
  float r,theta,phi;
  ARGST args;
  args.N=N;

  for(int i=0;i<N;i++){
    if( origDist(&(x[3*i])) > boundr ){
      outcount++;
      cnt=0;
      while(true){
	cnt++;
	r=boundr*mrand();
	theta=2*M_PI*mrand();
	phi=acos(2.0*mrand()-1.0);
	x[i*3]=r*sin(theta)*cos(phi);
	x[i*3+1]=r*sin(theta)*sin(phi);
	x[i*3+2]=r*cos(theta);
	if(LJpot(x,(void*)&args)-s->E < 10)
	  break;
	if(cnt>100){
	  printf("Error:unable to salt atom!\n");
	  exit(0);
	}
      }
    }
  }
  if(outcount>0)
    printf("salt cnt=%d\n",outcount);
}

