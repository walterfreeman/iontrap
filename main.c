#include "ionTrap.h"
#include<stdio.h>
#include <string>
#include <vector>
//==========================//
//           MAIN           //
//==========================//
int main(int argc, char *argv[]){
double Vac = atof(argv[1]);
double Vdc = atof(argv[2]);
double w = atof(argv[3]);
  ionTrap trap(Vac, Vdc, w);
    std::string rootFilename(argv[4]);
    trap.RenameRootFile(rootFilename);
    trap.BuildRootFile();
    //trap.PrintForAnime();
}
