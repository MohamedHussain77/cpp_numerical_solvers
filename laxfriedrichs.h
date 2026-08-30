#include "field.h"
#include "numericalmethod.h"
#ifndef LAXFRIEDRICHS_H
#define LAXFRIEDRICHS_H

class Laxfriedrichs1D: public NumericalMethod{
    private:
        double a;

    public:
        Laxfriedrichs1D(Field_u& u_param, double dx_param, double a_param):NumericalMethod(u_param, dx_param), a(a_param){}
        /*virtual void solve(double delt) override{
            for (size_t j{1}; j<(u.size()-1); ++j){
                u[j] = 0.5*(((1+(a*delt/delx))*u_old[j-1]) + ((1-(a*delt/delx))*(u_old[j+1])));
        }}*/
       virtual void apply(std::vector<double>& rhs, double dt) override{
            for (size_t j{1}; j<(u.size()-1); ++j){
                rhs[j] = ((0.5/dt)*(u[j+1] + u[j-1]) - (1.0/dt)*u[j])
                            - (0.5*a/dx)*(u[j+1] - u[j-1]);
        }            
        }
}; 
#endif