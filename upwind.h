#include "field.h"
#include "numericalmethod.h"
#ifndef UPWIND_H
#define UPWIND_H

class upwind: public NumericalMethod{

    public:
        upwind(Field_u& u_param, Field_u& u_old_param, double a_param, double delt_param, double delx_param):NumericalMethod(u_param,u_old_param, a_param, delt_param, delx_param){}
        virtual void solve() override{
            for (size_t j{1}; j<(u.size()-1); ++j){
                    u[j] = u_old[j] - (a*delt/delx)*(u_old[j] - u_old[j-1]);                   
                }
        }
}; 
#endif