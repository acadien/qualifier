#include "stdlib.h"
#include "stdio.h"
#include "gsl/gsl_math.h"
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

  //Prepare optimizer
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
    MCstep(&s,(void*)&args,true);
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

  //Test acceptance rate
  int accept=0;
  int testLoop=1000;
  for(int i;i<testLoop;i++){
    if(MCstep(&s,(void*)&args,false))
      accept++;
    basinPowell(&s,ftol,LJpotPunish,(void*)&args);
    printStateEnergy(&s,fp);
    printStateBounds(&s,fp);
  }
  printf("acceptance rate=%3.3f\n",(float)accept/testLoop);
  printStateBounds(&s,fp);
  
  FILE *logf;
  logf = fopen("finalstate.dat","w");
  printState(&s,logf);
  return 0;
}

bool MCstep(state* s, void* args,bool silent,float MCT, queue<int> accepts){
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
  float weight=exp(sprime.E/(float)s->N/MCT);
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
  return accept;
}

