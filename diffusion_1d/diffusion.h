#include "solver.h"
#ifndef DIFFUSION_H
#define DIFFUSION_H
class Diffusion : public Solver{ 
    private:
        double nu;
        double r;
    protected:
        virtual double getDeltaT() const override{
            return r * mesh.get_delx() * mesh.get_delx() / nu;
        }   


    public:
        Diffusion( const Mesh& mesh_param, Field& u_param, double nu_param, double r_param, double totaltime_param, NumericalMethod& method_param, TimeStepper& stepper_,
             BoundaryCondition& condition_param_l, BoundaryCondition& condition_param_r ) : Solver( mesh_param, u_param, totaltime_param, method_param, stepper_, condition_param_l, condition_param_r ), nu(nu_param), r(r_param) {}
            
};
#endif