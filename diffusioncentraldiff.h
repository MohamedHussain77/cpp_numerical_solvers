#include "field.h"
#include "numericalmethod.h"
#ifndef DIFFUSIONCENTRALDIFF_H
#define DIFFUSIONCENTRALDIFF_H

class Diffusioncentraldiff1D: public NumericalMethod{
    private:
        double a;

    public:
        Diffusioncentraldiff1D(Field_u& u_param,  double dx_param, double a_param):NumericalMethod(u_param, dx_param), a(a_param){}
        /*virtual void solve(double delt) override{
            for (size_t j{1}; j<(u.size()-1); ++j){
                u[j] = u_old[j] + (a*delt/(delx*delx))*(u_old[j+1]- 2 * u_old[j] + u_old[j-1]);
        }}*/
       virtual void apply(std::vector<double>& rhs, double dt) override{
            for (size_t j{1}; j<(u.size()-1); ++j){
                rhs[j] = + (a/(dx*dx))*(u[j+1]- 2 * u[j] + u[j-1]);
        }            
        }
}; 
#endif