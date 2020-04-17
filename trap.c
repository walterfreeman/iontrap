#include "iontrap-noroot.h"
#include<stdio.h>
#include <string>
#include <vector>
//==========================//
//           MAIN           //
//==========================//




int main(int argc, char *argv[]){
  double z0 = 3e-3;
  // Particle mass (in kg)
  double m = 2.28e-25;
  double q = -1.6e-19;

  if (argc < 8)
  {
    fprintf(stderr,"Usage: <this> <ac frequency in GHz> <az_min> <az_max> <az_step> <qz_min> <qz_max> <qz_step> <dt> <max cycles>\n");
    exit(1);
  }
  double w = atof(argv[1]); 
  double azmin = atof(argv[2]); 
  double azmax = atof(argv[3]); 
  double azstep= atof(argv[4]); 
  double qzmin = atof(argv[5]); 
  double qzmax = atof(argv[6]); 
  double qzstep= atof(argv[7]); 
  double dt = atof(argv[8]);
  double maxcycles = atof(argv[9]);
  int azcount=0;
  int qzcount=0;
  omp_set_num_threads(1);
  #pragma omp parallel 
  for (double az=azmin; az<azmax; az+=azstep)
  {
	  fprintf(stderr,"%f percent done.\n",(az-azmin)/(azmax-azmin)*100);
    for (double qz=qzmin; qz<qzmax; qz+=qzstep)
    {
      double Vdc = -az*m*z0*z0*w*w*1e18/(16*q);
      double Vac =  qz*m*z0*z0*w*w*1e18/(8*q);
      ionTrap trap(Vac, Vdc, w);
      trap.assess_stability(azcount, qzcount, dt, maxcycles);
      qzcount++;
    }
    azcount++; qzcount=0;
  }
}
