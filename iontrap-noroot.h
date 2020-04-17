#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include "vector.h"
#include <omp.h>



class ionTrap {
  public:

    double timer[100];

    void starttimer(int t)
    {
      timer[t]=omp_get_wtime();
    }

    // and this one stops it, returning how long it's been in microseconds since starttimer was called
    double checktimer(int t)
    {
      return omp_get_wtime()-timer[t];
    }

    double Vac;
    double Vdc;
    double w;
    double w_real = w*1e9;
    // For the Drag force (in kg/s)
    //double drag = 0;
    double drag = 3.47e-16;
    // For the Brownian kick (in J)
    double kT = 0;
    //double kT = 1.2e-21;
    // For the eletric force
    double q = -1.6e-19;
    double z0 = 3e-3;
    // Particle mass (in kg)
    double m = 2.28e-25;
    // Time step (in s)
    double dt = 1e-11;
    double t = 0;
    // Brownian kick
    double kickvar = 2*kT*drag/dt;
    double a_z = -16*Vdc*q/(m*z0*z0*w_real*w_real);
    double q_z = 8*Vac*q/(m*z0*z0*w_real*w_real);
    // Position and velocity (in m and m/s)
    vector r;
    vector v;
    vector Fext = vector(2e3*m, 1e3*m, 1e3*m); // external force
    double timeCounter = 0;
    double sumspeed;
    double speedFinal;
    double costhetasum=0;
    int costhetacount=0;
    std::string rootOutputFile = "trapOutput.root";
    //------------------------------------------//
    // Class methods
    ionTrap(double Vac, double Vdc, double w): Vac(Vac), Vdc(Vdc), w(w) {
      r=vector(0.000, 0.0000, 0.000)*0;
      v=vector(-0.0, 0.0, -0.0);
    }
    ~ionTrap(){}

    // Randomly generating the sign of the Brownian force term
    double kick_sign(){
      int random = rand() % 100 + 1;
      if(random % 2 == 0){
	return sqrt(kickvar);
      }
      else{
	return -sqrt(kickvar);
      }
    }

    // Calculating the force vector the particle is under
    vector F(vector r, vector v, double t){
      //Gravity force
      //The eletric term
      vector Fe = q * vector((Vdc - Vac * cos (w_real * t)) / (z0*z0) * ( 2 * r.x),
	  (Vdc - Vac * cos (w_real * t)) / (z0*z0) * ( 2 * r.y),
	  (Vdc - Vac * cos (w_real * t)) / (z0*z0) * (-1 * r.z));
      //The drag term
      vector Fd = -1 * v * drag;
      //The Brownian term
      vector Fb = vector((kick_sign() * kickvar), (kick_sign() * kickvar), (kick_sign() * kickvar));
      costhetacount++;
      costhetasum += (Fe+Fd+Fb) * v / (mag(v) * mag(Fe+Fd+Fb));
	


//      printf("\n\ntime %e\nposition:\t %e %e %e\nvelocity:\t %e %e %e\n",t,r.x,r.y,r.z,v.x,v.y,v.z);
//      printf("Fe       :\t%e %e %e\n",Fe.x,Fe.y,Fe.z);
//      printf("Fb       :\t%e %e %e\n",Fb.x,Fb.y,Fb.z);
//      printf("Fd       :\t%e %e %e\n",Fd.x,Fd.y,Fd.z);
 //     printf("Fext     :\t%e %e %e\n",Fext.x,Fext.y,Fext.z);
      return Fe + Fd + Fb + Fext;
    }
    void ExecuteLeapFrog() {
      r += v * dt/2;
      v += F(r,v,t)/m * dt/2;
      t += dt/2;
    }



    double assess_stability(int azcount, int qzcount, double timestep, double maxcycles){
      dt = timestep;
      printf("!Parameters:\n");
      printf("! Vac = %e \t Vdc = %e \t w = %e \t dt = %e \t drag = %e \t kT = %e\n",Vac,Vdc,w_real,dt,drag,kT);
      printf("! a_z = %e \t q_z = %e\n",a_z,q_z);
      printf("!Sanity checks: period %e | drag decay time %e | timestep %e \n",(1/w_real),m/drag,dt);
      double timescale_freq = 1/w_real;
      double timescale_terminal = m/drag;
      int steps=0;
      starttimer(0);
      starttimer(1);
      double stability=-1;
      while (true) {

	// make some decisions about whether we are stable or not

	if (drag == 0) // we'll have to decide stability based on things other than terminal velocity, because we don't have a terminal velocity here
	{
	  double rad = mag(r);
	  if (rad > z0)
	  {
	    printf("particle has escaped at time %e\n",t);
	    stability=0;
	    break;
	  }
	  if (t > 1e4/w_real)
	  {
	  	printf("Runtime = %f sec; completed %d steps, rate %f steps per second; particle at radius %e; radratio %e\n",checktimer(1),steps,steps/checktimer(1),mag(r),mag(r) / (t * mag(Fext)/drag));
		printf("declared stable after t = %e\n",t);
	        stability=1;
		break;
	  }
	}
	if (drag > 0)
	{
	  double rad = mag(r);
	 if (rad > z0)
          {
            printf("particle has escaped at time %e\n",t);
            stability=0;
            break;
          }
	  if (t > maxcycles / w_real)
	  {
		  double radratio = mag(r) / (t * mag(Fext)/drag);
		  stability = 1-radratio;
		  printf("Ending after %d cycles: radius ratio is %e, performance: %.2e steps/sec\n",(int)maxcycles,radratio,steps/checktimer(1));
		  if (radratio > 1) stability=0;
		  break;
	  }
	}
	for (int i=0; i<1e0; i++)
	{
	  ExecuteLeapFrog();
	  steps++;
	}

	if (checktimer(0) > 0.1)
	{
	  printf("Runtime = %f sec; copmleted %d steps, rate %f steps per second; particle at radius %e; radratio %e\n",checktimer(1),steps,steps/checktimer(1),mag(r), mag(r) / (t * mag(Fext)/drag));
	  starttimer(0);
	}

      }
      printf("RESULTS Vac %e Vdc %e w %e steps %d dt %e dragtime %e cycletime %e drag %e kT %e az %e qz %e avgcostheta %e finalrad %e termrad %e azcount %d qzcount %d stability %e\n",Vac,Vdc,w_real,steps,dt,timescale_terminal, timescale_freq, drag, kT, a_z, q_z, costhetasum/costhetacount, mag(r), t * mag(Fext)/drag,azcount, qzcount, stability);
    return stability;
    }
};
