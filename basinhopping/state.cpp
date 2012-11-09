#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "gsl/gsl_math.h"
//mine
#include "random.h"
#include "constants.h"
#include "state.h"
#include "structure.h"
#include "potential.h"

using namespace std;

void allocState(state* s,int N){
  s->N=N;
  s->x=(float*)malloc(sizeof(float)*3*N);  
  for(int i=0;i<3;i++)
    s->com[0]=0.;
}

void initState(state* s,int N){
  allocState(s,N);

  int failcount=0;

  float r,theta,phi;
  float* x=s->x; //test location
  bool found=false;

  ARGST args;
  args.N=N;

  failcount=0;
  found=false;
  while(!found){
    found=true;
    if(failcount++ > N*10000){
      printf("Too many trials on initialization, try a larger bounding sphere.\n");
      exit(0);
    }    
    for(int i=0;i<N;i++){
      //Select some random coordinates
      r=boundr*mrand();
      theta=2*M_PI*mrand();
      phi=acos(2.0*mrand()-1.0);
      x[i*3]=r*sin(theta)*cos(phi);
      x[i*3+1]=r*sin(theta)*sin(phi);
      x[i*3+2]=r*cos(theta);
    }

    if(LJpot(x,(void*)&args)>5E8)
      found=false;
  }
  
  return;
}

void printState(state* s,FILE *fp){
  printStateEnergy(s,fp);
  printStateVolume(s,fp);
  printStateBounds(s,fp);
  fprintf(fp,"COORDINATES\n");
  for(int i=0;i<s->N;i++){
    fprintf(fp,"% 5.5f % 5.5f % 5.5f\n",s->x[3*i],s->x[3*i+1],s->x[3*i+2]);
  }
}

void printStateEnergy(state* s, FILE *fp){
  fprintf(fp,"ENERGY\n");
  fprintf(fp,"% 6.6f\n",s->E);
}

void printStateBounds(state *s, FILE *fp){
  float mnx=1E10,mny=1E10,mnz=1E10,mxx=-1E10,mxy=-1E10,mxz=-1E10;
  float *xs=s->x;
  for(int i=0;i<s->N;i++){
    if(xs[3*i]<mnx)
      mnx=xs[3*i];
    if(xs[3*i]>mxx)
      mxx=xs[3*i];
    if(xs[3*i+1]<mny)
      mny=xs[3*i+1];
    if(xs[3*i+1]>mxy)
      mxy=xs[3*i+1];
    if(xs[3*i+2]<mnz)
      mnz=xs[3*i+2];
    if(xs[3*i+2]>mxz)
      mxz=xs[3*i+2];
  }
  fprintf(fp,"BOUNDS\n");
  fprintf(fp,"delx= %6.6f | ",mxx-mnx);
  fprintf(fp,"dely= %6.6f | ",mxy-mny);
  fprintf(fp,"delz= %6.6f\n",mxz-mnz);
}

void printStateVolume(state *s, FILE *fp){
  float mnx=1E10,mny=1E10,mnz=1E10,mxx=-1E10,mxy=-1E10,mxz=-1E10;
  float *xs=s->x;
  for(int i=0;i<s->N;i++){
    if(xs[3*i]<mnx)
      mnx=xs[3*i];
    if(xs[3*i]>mxx)
      mxx=xs[3*i];
    if(xs[3*i+1]<mny)
      mny=xs[3*i+1];
    if(xs[3*i+1]>mxy)
      mxy=xs[3*i+1];
    if(xs[3*i+2]<mnz)
      mnz=xs[3*i+2];
    if(xs[3*i+2]>mxz)
      mxz=xs[3*i+2];
  }
  fprintf(fp,"VOLUME\n");
  fprintf(fp,"%6.6f\n",(mxx-mnx)*(mxy-mny)*(mxz-mnz));
}

void copyState(state* a, state* b){ //copies state a on to state b
  if(b->N != a->N)
    allocState(b,a->N);
  b->E=a->E;
  for(int i=0;i<3*a->N;i++)
    b->x[i]=a->x[i];
}

