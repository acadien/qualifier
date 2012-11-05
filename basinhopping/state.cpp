#include "stdlib.h"
#include "stdio.h"
#include "gsl/gsl_math.h"
//mine
#include "random.h"
#include "constants.h"
#include "state.h"

using namespace std;
using namespace State;

void initState(state* s,int N){
  s->N=N;
  s->x=(float*)malloc(sizeof(float)*3*N);

  int failcount=0;

  float r,theta,phi,u,v;
  float tx,ty,tz; //test location
  bool found;
  for(int i=0;i<N;i++){
    found=false;
    while(!found){
      if(failcount++ > N*1000){
	printf("Too many trials on initialization, try a larger bounding sphere.\n");
	exit(0);
      }
      found=true;

      //Select some random coordinates
      r=boundr*mrand();
      u=mrand();
      v=mrand();
      theta=2*M_PI*u;
      phi=acos(2.0*v-1.0);
      tx=r*sin(theta)*cos(phi);
      ty=r*sin(theta)*sin(phi);
      tz=r*cos(theta);

      //Test coordinates      
      for(int j=0;j<i;j++){
	if( dist(&(s->x[i]),&(s->x[j])) < mindist ){
	  found=false;
	  break;
	}
      }
    }

    //Set the coordinates
    //printf("%d\n",i);
    s->x[3*i]=tx;
    s->x[3*i+1]=ty;
    s->x[3*i+2]=tz;
    
  }
  
  return;
}

void printState(state* s,FILE *fp){
  fprintf(fp,"% 6.6f\n",s->E);
  for(int i=0;i<s->N;i++){
    fprintf(fp,"% 5.5f % 5.5f % 5.5f\n",s->x[3*i],s->x[3*i+1],s->x[3*i+2]);
  }
}

float dist(float *a, float* b){
  statex=a[0]-b[0];
  statey=a[1]-b[1];
  statez=a[2]-b[2];
  return sqrt(statex*statex + statey*statey + statez*statez);
}
