#include <stdlib.h>
#include <stdio.h>
#include <queue>
#include <string>
#include "gsl/gsl_math.h"
//mine
#include "constants.h"
#include "random.h"
#include "state.h"
#include "localMin.h"
#include "structure.h"
#include "basinhop.h"

using namespace std;

int main(int argc, char **argv){
  state s;
  state sprime;
  state sideal;
  state* basins;

  //Settings
  int aLen=100;  //how big should the window average be
  int nAtom=38; //how many atoms
  //Basin
  float ftol=0.05; //set the tolerance on basin finding algo
  //MC
  int initLoop=1000, hopLoop=5000; //MC loop lengths
  float MCT=0.8;                  //Monte Carlo Temperature
  float MCalpha=0.30;             //Monte Carlo initial jump length

  //Initialize random numbers (mersenne twist)
  initrng();

  //Load up ideal cluster
  char buf[5];
  sprintf(buf,"%d",nAtom);
  string idealfilename=buf;
  idealfilename.append("_ideal.dat");
  FILE *idealfp;
  idealfp = fopen((char*)idealfilename.c_str(),"r");
  initState(&sideal,nAtom);
  loadIdeal(&sideal,idealfp);

  //Initialize cluster
  initState(&s,nAtom);
  initState(&sprime,nAtom);
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
  printf("Initial\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printf("\n");

  //Initial cluster will be horrible, run it through some MC to reduce total energy
  for(int i=0;i<initLoop;i++){
    MCstep(&s,&sprime,(void*)&args,ftol,MCT,MCalpha,&accepts,&acceptAvg,true);
    printStateVolume(&s,fp);
    printStateEnergy(&s,fp);
    printStateBounds(&s,fp);
  }
  printf("Reduced Initial.\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printStateBounds(&s,fp);
  printf("\n");
  
  /*
  //Obtain your first basin!
  basinPowell(&s,ftol,LJpot,(void*)&args);
  printf("Powell Reduced Initial.\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printStateBounds(&s,fp);
  printf("\n");
  */

  //Let the basin hopping begin
  float msdnow;
  for(int i=0;i<hopLoop;i++){

    MCstep(&s,&sprime,(void*)&args,ftol,MCT,MCalpha,&accepts,&acceptAvg,true);

    copyState(&s,&(basins[i]));

    if(i>0){
      msdnow=msd(&basins[i],&basins[i-1]);
      basins[i].msd=msdnow;
    }
    msdnow=msd(&basins[i],&sideal);
    basins[i].msdIdeal=msdnow;

    printf("****************************\n");
    printf("%d\n",i);
    if(i>0)
      printf("msd ideal=%f stepsize=%f acceptRatio=%f\n",msdnow,MCalpha,acceptAvg);
    printStateEnergy(&s,fp);
    printStateBounds(&s,fp);
    printf("****************************\n");
  }

  //Reoptimize with higher accuracy
  ftol=1e-4;
  float estart,eend;
  for(int i=0;i<hopLoop;i++){
    estart=basins[i].E;
    copyState(&(basins[i]),&sprime);
    basinPowell(&sprime,ftol,LJpot,(void*)&args);
    if(sprime.E<basins[i].E)
      printState(&sprime,logf);
    else
      printState(&(basins[i]),logf);
  }

  msdnow=msd(&basins[0],&basins[hopLoop-1]);
  printf("Total MSD=%f\n",msdnow);
  for(int i=0;i<hopLoop;i++)
    freeState(&(basins[i]));
  free(basins);
  freeState(&sprime);
  freeState(&sideal);
  freeState(&s);
  
  return 0;
}

void MCstep(state* s, state* sprime,void* args,float ftol, float& MCT, float& MCalpha, std::queue<int>* accepts, float* acceptAvg, bool silent){

  int cnt=0;
  float alphaStep=0.0002,alphaRatio=0.99;
  ARGST args2;
  args2.N=s->N;

  state sp=*sprime;

  while(true){
    cnt++;

    //Squeeze the structure to be more "cube-like" 
    cubify(s);

    //Salt atoms that are outside of sphere boundary
    salt(s);

    //Step out of local minimum!
    for(int i=0;i<sp.N;i++){
      sp.x[3*i+1] = s->x[3*i+1]+(mrand()-0.5)*2.0*MCalpha;
      sp.x[3*i+2] = s->x[3*i+2]+(mrand()-0.5)*2.0*MCalpha;
      sp.x[3*i+3] = s->x[3*i+3]+(mrand()-0.5)*2.0*MCalpha;
    }
    basinPowell(&sp,ftol,LJpot,(void*)&args2);
    sp.iters=0;
    //sp.E=LJpot(sp.x,args);

    //Calculate the Metropolis Criterion
    float weight=exp( -(sp.E - s->E) / MCT );
    if(!silent)
      printf("old:%4.4f new:%4.4f | expdelE=%4.4f\n",s->E,sp.E,weight);

    //Monte-Carlo action bam-pow
    bool accept=false;
    if(sp.E < s->E){
      //if(!silent)
      //  printf("lower energy\n");
      copyState(&sp,s);
      accept=true;
    }else 
      if(mrand()<weight){
	//if(!silent)
	//printf("higher energy\n");
	copyState(&sp,s);
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
    if(accept)  
      MCalpha/=alphaRatio;    
    else if (MCalpha > alphaStep)
      MCalpha*=alphaRatio;
    
    if(!silent){
      printf("inside: acceptAvg=%f MCT=%f           MCalpha=%f\n",*acceptAvg,MCT,MCalpha);
      //printf("===========================\n\n");
    }
    if(accept)
      break;
  }
  //  if(!silent)
  printf("cnt=%d\n",cnt);
}


void loadIdeal(state* sideal,FILE* ifile){
  int N=sideal->N,dummy;
  float x,y,z;
  for(int i=0;i<N;i++){
    dummy=fscanf(ifile,"%f %f %f\n",&x,&y,&z);
    sideal->x[3*i+1]=x;
    sideal->x[3*i+1]=y;
    sideal->x[3*i+1]=z;
  }
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
