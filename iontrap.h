#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include "vector.h"
#include "TFile.h"
#include "TTree.h"
class ionTrap {
public:
  double Vac;
  double Vdc;
  double w;
  double w_real = w*1e9;
  // For the Drag force (in kg/s)
  double drag = 0;
  //double drag = 3.47e-16;
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
  double timeCounter = 0;
  double sumspeed;
  double speedFinal;
 std::string rootOutputFile = "trapOutput.root";
  //------------------------------------------//
  // Class methods
  ionTrap(double Vac, double Vdc, double w): Vac(Vac), Vdc(Vdc), w(w) {
    r.reassign(0.001, 0.0005, 0.001);
    v.reassign(-0.5, 0.5, -0.1);
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
    vector Fg = vector(0, -9.8*m, 0);
    //The eletric term
    vector Fe = q * vector((Vdc - Vac * cos (w_real * t)) / (z0*z0) * ( 2 * r.x),
                (Vdc - Vac * cos (w_real * t)) / (z0*z0) * ( 2 * r.y),
                (Vdc - Vac * cos (w_real * t)) / (z0*z0) * (-1 * r.z));
    //The drag term
    vector Fd = -1 * v * drag;
    //The Brownian term
    vector Fb = vector((kick_sign() * kickvar), (kick_sign() * kickvar), (kick_sign() * kickvar));
    return Fe + Fd + Fb + Fg;
  }
  void ExecuteLeapFrog() {
    r += v * dt/2;
    v += F(r,v,t)/m * dt/2;
    t += dt/2;
}
  void BuildRootFile(){
    // Defining vectors to be used in the root file
    std::vector<double> vector_X;
    std::vector<double> vector_Y;
    std::vector<double> vector_Z;
    std::vector<double> vector_vX;
    std::vector<double> vector_vY;
    std::vector<double> vector_vZ;
    std::vector<double> vector_speedFinal;
    std::vector<double> vector_Vac;
    std::vector<double> vector_Vdc;
    std::vector<double> vector_frequency;
    std::vector<double> vector_normR;
    std::vector<double> vector_normZ;
    // Creating the root file
    TFile rootOutput((const char*) rootOutputFile.c_str(), "RECREATE");
TTree inputTree("inputs", "inputs");
    inputTree.Branch("brownianKick", &kT, "brownianKick/D");
    inputTree.Branch("drag", &drag, "drag/D");
    inputTree.Branch("q", &q, "q/D");
    inputTree.Branch("z0", &z0, "z0/D");
    inputTree.Branch("m", &m, "m/D");
    inputTree.Branch("dt", &dt, "dt/D");
    TTree trapTree("trap", "trap");
    trapTree.Branch("X", &vector_X);
    trapTree.Branch("Y", &vector_Y);
    trapTree.Branch("Z", &vector_Z);
    trapTree.Branch("norm_R", &vector_normR);
    trapTree.Branch("norm_Z", &vector_normZ);
    trapTree.Branch("vX", &vector_vX);
    trapTree.Branch("vY", &vector_vY);
    trapTree.Branch("vZ", &vector_vZ);
    trapTree.Branch("Vac", &vector_Vac);
    trapTree.Branch("Vdc", &vector_Vdc);
    trapTree.Branch("w", &vector_frequency);
    trapTree.Branch("Final_speed", &vector_speedFinal);
 // Calculating kickvar
    kickvar = 2 * kT * drag/dt;
    while(t < (10000/w_real)){
      ExecuteLeapFrog();
      double sumspeed =+ norm(v);
    }
    double normR = sqrt(r.x*r.x + r.y*r.y);
    double normZ = sqrt(r.z*r.z);
    speedFinal = sumspeed/t;
    vector_X.push_back(r.x);
    vector_Y.push_back(r.y);
    vector_Z.push_back(r.z);
    vector_vX.push_back(v.x);
    vector_vY.push_back(v.y);
    vector_vZ.push_back(v.z);
    vector_speedFinal.push_back(speedFinal);
    vector_Vac.push_back(Vac);
    vector_Vdc.push_back(Vdc);
    vector_frequency.push_back(w_real);
    vector_normR.push_back(normR);
    vector_normZ.push_back(normZ);
  r.reassign(0.001, 0.0005, 0.001);
    v.reassign(-0.5, 0.5, -0.1);
    t = 0;
    // For the root file
    inputTree.Fill();
    inputTree.Write();
    trapTree.Fill();
    trapTree.Write();
    rootOutput.Close();
  }
  void PrintForAnime(){
    double time_factor = 1e-5; // fraction of realtime to run at
    int frameskip = 1./(60 * dt) * time_factor;
    printf("!Beginning simulation. Frameskip = %d\n",frameskip);
    printf("!Parameters:\n");
    printf("! Vac = %e \t Vdc = %e \t w = %e \t dt = %e \t drag = %e \t kT = %e\n",Vac,Vdc,w_real,dt,drag,kT);
    printf("! a_z = %e \t q_z = %e\n",a_z,q_z);
    printf("!Sanity checks: dimensionless period %e | drag decay constant %e\n",dt*w_real,drag/m * dt);
    int steps;
    while (true) {
      ExecuteLeapFrog();
 steps++;
      if (steps % frameskip == 0)
      {
        printf("ct3 0 %e %e %e %e\n", r.x, r.y, r.z, 1e-4);
        printf("F\n");
      }
    }
  }
  void RenameRootFile(std::string &rootFileName){
    this->rootOutputFile = rootFileName;
  }
};
