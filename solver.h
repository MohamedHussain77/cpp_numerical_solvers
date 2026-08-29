#include<vector>
#include "mesh.h"
#include "field.h"
#include "numericalmethod.h"
#include "boundarycondition.h"
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
        
        NumericalMethod& method;
        BoundaryCondition& condition_left;
        BoundaryCondition& condition_right;
        

    public:
        Solver_advection(const Mesh& mesh_param,Field_u& u_param, Field_u& u_old_param, double a_param, double cfl_param, 
            double totaltime_param, NumericalMethod& method_param, BoundaryCondition& condition_param_l, BoundaryCondition& condition_param_r): mesh(mesh_param), a(a_param), cfl(cfl_param), delt(cfl*mesh.get_delx()/std::abs(a)),
            totaltime(totaltime_param), numsteps(static_cast<int>(totaltime / delt)), u(u_param), u_old(u_old_param), method(method_param), condition_left(condition_param_l), condition_right(condition_param_r){}
        
        void advection(){
            u_old.copy_from(u);
            std::ofstream output("solution_ADVECTION_tophat.csv");
            double time = 0.0;
            output << "t,x,u\n";
            for (size_t j{0}; j<mesh.get_numcells(); ++j){
                    output << time << "," << mesh[j] << "," << u[j] << "\n";
                }
            for (size_t i{0}; i<numsteps; ++i){
                //double time = delt * i;
                method.solve();                    
                
                condition_left.apply();
                condition_right.apply();
                time += delt;
                for (size_t j{0}; j<mesh.get_numcells(); ++j){
                        output << time << "," << mesh[j] << "," << u[j] << "\n";
                    }
                //u[0]=0;
                //u[mesh.get_numcells()-1]=0;
                u_old.swap_btw(u);    
            }
            output.close();
        }
        int get_numsteps() const{
            return numsteps;
        }


};
#endif



