#include "field.h"
#include "numericalmethod.h"
#ifndef UPWIND_H
#define UPWIND_H

class upwind1D: public NumericalMethod{

    private:
        double a;

    public:
        upwind1D(double dx_param, double a_param):NumericalMethod(dx_param), a(a_param){}
        /*virtual void solve(double delt) override{
            for (size_t j{1}; j<(u.size()-1); ++j){
                    //u[j] = u_old[j] - (a*delt/delx)*(u_old[j] - u_old[j-1]);                   
                }
        }*/
        virtual void apply(Field& u, std::vector<double>& rhs, double dt) override{

            for (size_t j{1}; j< (u.size()-1); ++j){
                rhs[j] = - a * (u[j]-u[j-1]) / dx;
            }
        }
}; 
#endif