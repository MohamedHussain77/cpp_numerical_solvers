#include <iostream>
#include <vector>
#include <fstream>
#include "mesh.h"
#include "field.h"
#include "solver.h"

int main(){
    double L = 1;
    int numcells = 100;
    Mesh mesh(numcells, L);

    double cfl = 0.8;
    double a = 2.1;
    double totaltime = 0.5;
    
    Field_u u(mesh);
    u.initialize_u_tophat();
    Field_u u_old(mesh);

    Solver_advection solve(mesh, u, u_old, a, cfl, totaltime);

    std::cout << "del x: " << mesh.get_delx() << std::endl;
    std::cout << "numsteps: " << solve.get_numsteps() << std::endl;

    solve.solving_advection_Laxwendroff();
    return 0; 
}