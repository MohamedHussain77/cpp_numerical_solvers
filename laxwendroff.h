#include "field.h"
#include "numericalmethod.h"
#ifndef LAXWENDROFF_H
#define LAXWENDROFF_H

class Laxwendroff: public NumericalMethod{

    public:
        Laxwendroff(Field_u& u_param, Field_u& u_old_param, double a_param, double delt_param, double delx_param):NumericalMethod(u_param,u_old_param, a_param, delt_param, delx_param){}
        virtual void solve() override{
            for (size_t j{1}; j<(u.size()-1); ++j){
                u[j] = u_old[j] - ((0.5*(a*delt/delx))*(u_old[j+1] - u_old[j-1]))+ ((0.5*(a*delt/delx)*(a*delt/delx))*(u_old[j+1] - 2*u_old[j] + u_old[j-1]));
        }}
}; 
#endif
