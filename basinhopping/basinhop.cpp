#include "stdlib.h"
#include "stdio.h"
#include <queue>
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
  state s;
  state* basins;

  //Settings
  int aLen=500;  //how big should the window average be
  int nAtom=38; //how many atoms
  //Basin
  float ftol=0.2; //set the tolerance on basin finding algo
  //MC
  int initLoop=100, hopLoop=1000; //MC loop lengths
  float MCT=0.8;                  //Monte Carlo Temperature
  float MCalpha=0.36;             //Monte Carlo initial jump length

  //Initialize random numbers (mersenne twist)
  initrng();

  //Initialize cluster
  initState(&s,nAtom);
  basins=(state*)malloc(sizeof(state)*hopLoop);
  for(int i=0;i<hopLoop;i++)
    allocState(&(basins[i]),nAtom);
  ARGST args;
  args.N=s.N;

  //Initalize logging
  s.E=LJpotPunish(s.x,(void*)&args);
  FILE* fp=stdout;
  FILE *logf;
  logf = fopen("finalstate.dat","w");

  //Initialize acceptance queue for dynamically altering variation parameter
  std::queue<int> accepts;
  float acceptAvg=0.5;
  for(int i=0;i<aLen/2;i++){
    accepts.push(1);
    accepts.push(0);
  }

  //Prepare optimizers
  initPowell(s.N*3);
  printf("Initial\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printf("\n");

  /*
  //Initial cluster will be horrible, run it through some MC to reduce total energy
  for(int i=0;i<initLoop;i++)
    MCstep(&s,(void*)&args,MCT,MCalpha,&accepts,&acceptAvg);
  printf("Reduced Initial.\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printStateBounds(&s,fp);
  printf("\n");
  */

  //Obtain your first basin!
  basinPowell(&s,ftol,LJpot,(void*)&args);
  printf("Powell Reduced Initial.\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printStateBounds(&s,fp);
  printf("\n");

  //Let the basin hopping begin
  //  MCalpha=0.1;             //Monte Carlo initial jump length
  float msdnow;
  for(int i=0;i<hopLoop;i++){

    MCstep(&s,(void*)&args,ftol,MCT,MCalpha,&accepts,&acceptAvg,true);

    //basinPowell(&s,ftol,LJpot,(void*)&args);

    copyState(&s,&(basins[i]));

    if(i>0)
      msdnow=msd(&basins[i],&basins[i-1]);

    fprintf(logf,"NAtoms\n%d\n",nAtom);
    printState(&s,logf);

    fprintf(logf,"MSD\n");
    if(i==0)
      fprintf(logf,"0\n");
    else
      fprintf(logf,"%f\n",msdnow);

    printf("****************************\n");
    printf("%d\n",i);
    if(i>0)
      printf("msd=%f stepsize=%f acceptRatio=%f\n",msdnow,MCalpha,acceptAvg);
    printStateEnergy(&s,fp);
    printStateBounds(&s,fp);
    printf("****************************\n");

  }

  printStateBounds(&s,fp);
  
  return 0;
}

void MCstep(state* s, void* args,float ftol, float& MCT, float& MCalpha, std::queue<int>* accepts, float* acceptAvg, bool silent){

  int cnt=0;
  state sprime;
  float alphaStep=0.0002,alphaRatio=0.99;
  initState(&sprime,s->N);
  initPowell(s->N*3);

  while(true){
    cnt++;
    //Salt atoms that are outside of sphere boundary
    salt(s);

    //Step out of local minimum!
    for(int i=0;i<sprime.N;i++){
      sprime.x[3*i] = s->x[3*i]+(mrand()-0.5)*2.0*MCalpha*0.5;
      sprime.x[3*i+1] = s->x[3*i+1]+(mrand()-0.5)*2.0*MCalpha*0.5;
      sprime.x[3*i+2] = s->x[3*i+2]+(mrand()-0.5)*2.0*MCalpha*0.5;
      /*
      sprime.x[3*i] = s->x[3*i]+mnormrand(MCalpha/10.);
      sprime.x[3*i+1] = s->x[3*i+1]+mnormrand(MCalpha/10.);
      sprime.x[3*i+2] = s->x[3*i+2]+mnormrand(MCalpha/10.);
      */
    }
    
    basinPowell(&sprime,ftol,LJpot,args);

    //sprime.E=LJpot(sprime.x,args);

    float weight=exp( -(sprime.E - s->E) / MCT );
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
      if(mrand()<weight){
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
    if(*acceptAvg > 0.52)
      MCalpha/=alphaRatio;
    else if (MCalpha > alphaStep && *acceptAvg < 0.48)
      MCalpha*=alphaRatio;
    
    if(!silent){
      printf("inside: acceptAvg=%f MCT=%f           MCalpha=%f\n",*acceptAvg,MCT,MCalpha);
      printf("===========================\n\n");
    }
    if(accept)
      break;
  }
  //  if(!silent)
  printf("cnt=%d\n",cnt);
}

void resetWindow(std::queue<int>* accepts,float* acceptAvg, int aLen){
  for(int i=0;i<aLen/2;i++){
    accepts->pop();
    accepts->pop();
    accepts->push(1);
    accepts->push(0);
  }
  *acceptAvg=0.5;
}
