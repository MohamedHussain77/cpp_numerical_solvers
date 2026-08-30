#include "solver.h"
#ifndef ADVECTION_H
#define ADVECTION_H
class Advection : public Solver{ 
    private:
        double a;
        double cfl;
    protected:
        virtual double getDeltaT() const override{
            return cfl * mesh.get_delx() / std::abs(a);
        }   


    public:
        Advection( const Mesh& mesh_param, Field_u& u_param, double a_param, double cfl_param, double totaltime_param, NumericalMethod& method_param, TimeStepper& stepper_,
             BoundaryCondition& condition_param_l, BoundaryCondition& condition_param_r ) : Solver(mesh_param, u_param, totaltime_param, method_param, stepper_, condition_param_l, condition_param_r ), a(a_param), cfl(cfl_param) {}
            
};
#endif