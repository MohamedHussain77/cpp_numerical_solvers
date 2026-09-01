#ifndef SOLVER_H
#define SOLVER_H
#include<vector>
#include "mesh.h"
#include "field.h"
#include "eulerflux1d.h"
#include "numericalmethod.h"
#include "boundarycondition.h"
#include "timestepper.h"
#include "numericalflux.h"
#include<fstream>


class Solver{
    protected:
        const Mesh& mesh;
        Field& U;
        //Field& u_old;
        //double a;
        double cfl;
        double gamma;
        //double delt;
        double totaltime;
        //int numsteps; 
        NumericalFlux& Flux;
        NumericalMethod& method;
        TimeStepper& stepper;
        BoundaryCondition& condition_left;
        BoundaryCondition& condition_right;

        

    public:
        
        Solver(const Mesh& mesh_param,Field& u_param, double cfl_, double gamma_,
            double totaltime_param, NumericalFlux& Flux_, NumericalMethod& method_param, TimeStepper& stepper_, BoundaryCondition& condition_param_l, BoundaryCondition& condition_param_r): mesh(mesh_param),
            U(u_param), cfl(cfl_), gamma(gamma_), totaltime(totaltime_param), Flux(Flux_), method(method_param), stepper(stepper_),condition_left(condition_param_l), condition_right(condition_param_r){}
        
        virtual ~Solver() = default;

        virtual double getDeltaT(const Field& U) const {
            double Smax = 0;
            for (std::size_t i{0}; i < U.size(); ++i){
                double S = std::abs(U[i].velocity()) + U[i].soundspeed(gamma);
                if (S > Smax){
                    Smax = S;
                }
            }
            return cfl * mesh.get_delx() / Smax;
        }

        virtual void timeIntegrate(){
            std::vector<ConservedState1D> rhs(mesh.get_numcells());
            std::vector<EulerFlux1D> F(mesh.get_numcells()+1);
            //u_old.copy_from(u);
            //double delt = getDeltaT(U, gamma);
            std::ofstream output("solution.csv");
            double time = 0.0;
            output << "t,x,rho,u,p\n";
            for (size_t j{0}; j<mesh.get_numcells(); ++j){
                    output << time << "," << mesh[j] << "," << U[j].rho << "," << U[j].velocity() << "," << U[j].pressure(gamma) << "\n";
                }
            while (time < totaltime){
                double dt = getDeltaT(U);
                dt = std::min(dt, totaltime - time);
                method.computeRHS(U, F, rhs, dt);
                stepper.step(U, rhs, dt);


               
                //method.apply(u, rhs, dt);
                //stepper.step(u, rhs, dt);    
                               
                
                //condition_left.apply();
                //condition_right.apply();
                time += dt;
                for (size_t j{0}; j<mesh.get_numcells(); ++j){
                    output << time << "," << mesh[j] << "," << U[j].rho << "," << U[j].velocity() << "," << U[j].pressure(gamma) << "\n";
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
