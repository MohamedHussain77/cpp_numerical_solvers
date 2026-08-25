#include<vector>
#include "mesh.h"
#include "field.h"
#include<fstream>
#ifndef SOLVER_H
#define SOLVER_H

class Solver_advection{
    private:
        const Mesh& mesh;
        Field_u& u;
        Field_u& u_old;
        double a;
        double cfl;
        double delt;
        double totaltime;
        int numsteps;        

    public:
        Solver_advection(const Mesh& mesh_param,Field_u& u_param, Field_u& u_old_param, double a_param, double cfl_param, 
            double totaltime_param): mesh(mesh_param), a(a_param), cfl(cfl_param), delt(cfl*mesh.get_delx()/a),
            totaltime(totaltime_param), numsteps(totaltime/delt), u(u_param), u_old(u_old_param){}
        
        void solving_advection_firstorderupwind(){
            u_old.copy_from(u);
            std::ofstream output("solution1.csv");
            output << "t,x,u\n";
            for (size_t i{0}; i<numsteps; ++i){
                double time = delt * i;
                for (size_t j{0}; j<mesh.get_numcells(); ++j){
                    output << time << "," << mesh[j] << "," << u_old[j] << "\n";
                }
                for (size_t j{1}; j<(mesh.get_numcells()-1); ++j){
                    u[j] = u_old[j] - (cfl)*(u_old[j] - u_old[j-1]);                    
                }
            u[0]=0;
            u[mesh.get_numcells()-1]=0;
            u_old.swap_btw(u);    
            }
            output.close();
        }
        int get_numsteps() const{
            return numsteps;
        }


};
#endif




