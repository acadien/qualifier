#include "stdlib.h"
#include "stdio.h"
#include "gsl/gsl_math.h"
//mine
#include "constants.h"
#include "random.h"
#include "potential.h"
#include "state.h"
#include "localMin.h"
#include "basinhop.h"

using namespace std;

int main(int argc, char **argv){
  state s;

  //Initialize RNG
  initrng();

  //Initialize cluster
  initState(&s,100);
  ARGST args;
  args.N=s.N;
  FILE* fp=stdout;

  //Initialize optimizer
  initPowell(s.N*3);

  //Start with some monte carlo stuff...
  
  /*
  //TEST: check energy of cluster  
  s.E=LJpot(s.x,(void*)&args);
  printState(&s,fp);

  //Optimize and print

  basinPowell(&s,LJpot,(void*)&args);
  //  s.E=LJpot(s.x,(void*)&args);
  printState(&s,fp);
  */
  return 0;
}


