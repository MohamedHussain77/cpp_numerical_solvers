#include "field.h"
#include "numericalmethod.h"
#ifndef LAXWENDROFF_H
#define LAXWENDROFF_H

class Laxwendroff1D: public NumericalMethod{
    private:
        double a;

    public:
        Laxwendroff1D(double dx_param, double a_param):NumericalMethod(dx_param), a(a_param){}
        /*virtual void solve(double delt) override{
            for (size_t j{1}; j<(u.size()-1); ++j){
                u[j] = u_old[j] - ((0.5*(a*delt/delx))*(u_old[j+1] - u_old[j-1]))+ ((0.5*(a*delt/delx)*(a*delt/delx))*(u_old[j+1] - 2*u_old[j] + u_old[j-1]));
        }}*/
        virtual void apply(Field& u, std::vector<double>& rhs, double dt) override{
            for (size_t j{1}; j<(u.size()-1); ++j){
                rhs[j] = - ((0.5*(a/dx))*(u[j+1] - u[j-1]))+ ((0.5*(a*dt/dx)*(a/dx))*(u[j+1] - 2*u[j] + u[j-1]));
        }            
        }
}; 
#endif
