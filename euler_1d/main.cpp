#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <numbers>
#include "mesh.h"
#include "field.h"
#include "uniforminitialcondition.h"
#include "sodshocktube.h"
#include "boundarycondition.h"
#include "numericalflux.h"
#include "numericalmethod.h"
#include "solver.h"
#include "rk2stepper.h"
#include "roeflux.h"
#include "laxfriedrichsflux.h"
#include "hllcflux.h"
#include "reconstruction.h"
#include "limiter.h"



int main(){
    double L = 1;
    int numcells = 250;
    double totaltime = 0.2;
    double cfl_max = 0.5;

    Mesh mesh(numcells, L);

    Field U(mesh);
    //Initialize internal field
    //double rho_internal = 1;
    //double u_internal = 0;
    //double p_internal = 1;
    double gamma = 1.4;
    //UniformInitialCondition initial(U, mesh, rho_internal, u_internal, p_internal, gamma);
    SodShockTubeTest initial(U, mesh, 1, 0, 1, 0.5, 0.125, 0, 0.1, gamma);
    initial.initialize();

    zeroGradient left(BoundarySide::Left);
    zeroGradient right(BoundarySide::Right);

    //LaxFriedrichsFlux flux(mesh.get_delx(),gamma);
    //RoeFlux flux(gamma);
    HLLCFlux flux(gamma);
    VanLeer vanleer;
    Reconstruction reconstruct(vanleer);
    NumericalMethod method(mesh, flux, left, right, reconstruct);
    RK2Stepper stepper(method);
    std::string filename("solution_hllc_vanleer.csv");

    Solver Euler(mesh, U, cfl_max, gamma, totaltime, flux, method, stepper, left, right, filename);
    Euler.timeIntegrate();



    


    return 0;

}
