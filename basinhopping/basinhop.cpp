#include "stdlib.h"
#include "stdio.h"
#include "gsl/gsl_math.h"
#include <queue>
//mine
#include "constants.h"
#include "random.h"
#include "potential.h"
#include "state.h"
#include "localMin.h"
#include "structure.h"
#include "basinhop.h"

using namespace std;

int main(int argc, char **argv){
  float ftol;
  state s;

  //Initialize RNG
  initrng();

  //Initialize cluster
  initState(&s,35);
  ARGST args;
  args.N=s.N;
  s.E=LJpot(s.x,(void*)&args);
  FILE* fp=stdout;

  //Initialize acceptance queue for dynamically altering temperature
  std::queue<int> accepts;
  int aLen=50;//how big should the window average be
  float acceptAvg=0.5;
  for(int i=0;i<aLen/2;i++){
    accepts.push(1);
    accepts.push(0);
  }

  //Prepare optimizer
  float MCT=0.8;
  initPowell(s.N*3);
  printf("Initial\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printf("\n");

  /*
  //TEST - must bound powell's method
  ftol=1E1;
  basinPowell(&s,ftol,LJpotPunish,(void*)&args);
  printf("Powell Reduced Initial.\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printStateBounds(&s,fp);
  printf("\n");
  */
    
  //Initial cluster will be horrible, run it through some MC to reduce total energy
  int initLoop=1000;
  for(int i;i<initLoop;i++){
    //    s.E=LJpot(s.x,(void*)&args);
    MCstep(&s,(void*)&args,&MCT,(std::queue<int>*)&accepts,&acceptAvg);
  }
  printf("Reduced Initial.\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printStateBounds(&s,fp);
  printf("\n");

  //Obtain your first basin!
  ftol=1.0;
  basinPowell(&s,ftol,LJpotPunish,(void*)&args);
  printf("Powell Reduced Initial.\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printStateBounds(&s,fp);
  printf("\n");

  //Let the basin hopping begin
  int testLoop=1000;
  for(int i;i<testLoop;i++){
    MCstep(&s,(void*)&args,&MCT,&accepts,&acceptAvg);
    printf("acceptance rate=%3.3f,  temperature=%4.4f\n",acceptAvg,MCT);
    basinPowell(&s,ftol,LJpotPunish,(void*)&args);
    printStateEnergy(&s,fp);
    printStateBounds(&s,fp);
  }

  printStateBounds(&s,fp);
  
  FILE *logf;
  logf = fopen("finalstate.dat","w");
  printState(&s,logf);
  return 0;
}

void MCstep(state* s, void* args,float* MCT, std::queue<int>* accepts, float* acceptAvg){
  bool silent=true;
  float origin[3] = {0., 0., 0.};
  state sprime;
  allocState(&sprime,s->N);


  //  if(sprime.x[10]!=s->x[10] || sprime.E!=s->E){
  //    printf("****ERROR COPY**** %f %f\n",sprime.E,s->E);
  //    exit(0);
  //  }
  while(true){
    copyState(s,&sprime); //copy s to sprime

    //Step out of local minimum!
    for(int i;i<sprime.N;i++){
      //    while(true){
      sprime.x[3*i] += (mrand()-0.5)*2.0*MCalpha*0.5;
      sprime.x[3*i+1] += (mrand()-0.5)*2.0*MCalpha*0.5;
      sprime.x[3*i+2] += (mrand()-0.5)*2.0*MCalpha*0.5;
      //if(dist(origin,
    }
    sprime.E=LJpotPunish(sprime.x,args);

    //  float weight=exp(-(s->E-sprime.E)/MCT/s->N);
    float weight=exp(sprime.E / (float)s->N / *MCT);
    if(!silent)
      printf("old:%4.4f new:%4.4f | expdelE=%4.4f\n",s->E,sprime.E,weight);

    //Monte-Carlo action bam-pow
    bool accept=false;
    if(sprime.E < s->E){
      if(!silent)
	printf("lower energy\n");
      copyState(&sprime,s);
      accept=true;
    }else 
      if(weight<FLT_MAX && mrand()>=weight){
	if(!silent)
	  printf("higher energy\n");
	copyState(&sprime,s);
	accept=true;
      }else{
	if(!silent)
	  printf("higher energy didn't take\n");
      }

    //Update the acceptance queue and average acceptance
    *acceptAvg-=accepts->front()/(float)accepts->size();
    accepts->pop();
    if(accept){
      accepts->push(1);
      *acceptAvg+=1.0/(float)accepts->size();
    }
    else
      accepts->push(0);

    //Update the Monte-Carlo Temperature according to acceptance
    if(*acceptAvg>0.52)
      *MCT+=0.1;
    else if (*acceptAvg < 0.48)
      *MCT-=0.1;
  }
}

