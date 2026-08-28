#include <iostream>
#include <vector>
#include <fstream>
#include "mesh.h"
#include "field.h"
#include "solver.h"
#include "numericalmethod.h"
#include "upwind.h"
#include "laxwendroff.h"
#include "dirichlet.h"
#include "neumann.h"
#include "robin.h"

int main(){
    double L = 1.5;
    int numcells = 150;
    Mesh mesh(numcells, L);

    double cfl = 0.4;
    double a = 2.5;
    double totaltime = 0.5;
    
    Field_u u(mesh);
    u.initialize_u_tophat();
    Field_u u_old(mesh);
    Laxwendroff laxwendroff (u, u_old);
    Dirichlet left(u, 0, BoundarySide::Left, mesh.get_delx());
    Robin right(u, 0, BoundarySide::Right, mesh.get_delx(), 0.2, 0.8); //alpha * u + beta du/dx = C 


    Solver_advection solve(mesh, u, u_old, a, cfl, totaltime, laxwendroff, left, right);

    std::cout << "del x: " << mesh.get_delx() << std::endl;
    std::cout << "numsteps: " << solve.get_numsteps() << std::endl;

    solve.advection();
    return 0; 
}