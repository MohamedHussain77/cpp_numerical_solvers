#include "field.h"
#include "numericalmethod.h"
#ifndef UPWIND_H
#define UPWIND_H

class upwind: public NumericalMethod{

    public:
        upwind(Field_u& u_param, Field_u& u_old_param):NumericalMethod(u_param,u_old_param){}
        virtual void solve(size_t j, double cfl) override{
           u[j] = u_old[j] - (cfl)*(u_old[j] - u_old[j-1]);
        }
}; 
#endif