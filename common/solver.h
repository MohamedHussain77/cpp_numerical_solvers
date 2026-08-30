#include<vector>
#include "mesh.h"
#include "field.h"
#include "numericalmethod.h"
#include "boundarycondition.h"
#include "timestepper.h"
#include<fstream>
#ifndef SOLVER_H
#define SOLVER_H

class Solver{
    protected:
        const Mesh& mesh;
        Field& u;
        //Field& u_old;
        //double a;
        //double cfl;
        //double delt;
        double totaltime;
        //int numsteps; 
        
        NumericalMethod& method;
        TimeStepper& stepper;
        BoundaryCondition& condition_left;
        BoundaryCondition& condition_right;

        virtual double getDeltaT() const = 0;

    public:
        
        Solver(const Mesh& mesh_param,Field& u_param, //Field& u_old_param, 
            double totaltime_param, NumericalMethod& method_param, TimeStepper& stepper_, BoundaryCondition& condition_param_l, BoundaryCondition& condition_param_r): mesh(mesh_param),
            u(u_param), totaltime(totaltime_param), method(method_param), stepper(stepper_),condition_left(condition_param_l), condition_right(condition_param_r){}
        
        virtual ~Solver() = default;

        virtual void timeIntegrate(){
            std::vector<double> rhs(mesh.get_numcells());
            //u_old.copy_from(u);
            double delt = getDeltaT();
            std::ofstream output("solution.csv");
            double time = 0.0;
            output << "t,x,u\n";
            for (size_t j{0}; j<mesh.get_numcells(); ++j){
                    output << time << "," << mesh[j] << "," << u[j] << "\n";
                }
            while (time < totaltime){
                double dt = std::min(delt, totaltime - time);
               
                method.apply(u, rhs, dt);
                stepper.step(u, rhs, dt);    
                               
                
                condition_left.apply();
                condition_right.apply();
                time += dt;
                for (size_t j{0}; j<mesh.get_numcells(); ++j){
                        output << time << "," << mesh[j] << "," << u[j] << "\n";
                    }
                //u[0]=0;
                //u[mesh.get_numcells()-1]=0;
                //u_old.swap_btw(u);    
            }
            output.close();
        }
        /*int get_numsteps() const{
            return numsteps;
        
        }*/
};
#endif
