#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <numbers>
#include "mesh.h"
#include "field.h"
#include "solver.h"
#include "numericalmethod.h"
#include "upwind.h"
#include "laxwendroff.h"
#include "dirichlet.h"
#include "neumann.h"
#include "robin.h"
#include "initialcondition.h"
#include "tophat.h"
#include "sinewave.h"
#include "gaussian.h"

int main(){
    double L = 1;
    int numcells = 150;
    Mesh mesh(numcells, L);

    double cfl = 0.5;
    double a = 3;
    double totaltime = 0.5;
    
    Field_u u(mesh);
    SineWave sine(u, mesh, 0.8, 4*std::numbers::pi,0,0);
    sine.initialize();
    //u.initialize_u_tophat();
    Field_u u_old(mesh);
    Laxwendroff laxwendroff (u, u_old);
    Dirichlet left(u, 0, BoundarySide::Left, mesh.get_delx());
    Neumann right(u, 0, BoundarySide::Right, mesh.get_delx()); //alpha * u + beta du/dx = C 


    Solver_advection solve(mesh, u, u_old, a, cfl, totaltime, laxwendroff, left, right);

    std::cout << "del x: " << mesh.get_delx() << std::endl;
    std::cout << "numsteps: " << solve.get_numsteps() << std::endl;

    solve.advection();
    return 0; 
}