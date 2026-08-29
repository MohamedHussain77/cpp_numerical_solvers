#include "field.h"
#include "numericalmethod.h"
#ifndef LAXFRIEDRICHS_H
#define LAXFRIEDRICHS_H

class Laxfriedrichs: public NumericalMethod{

    public:
        Laxfriedrichs(Field_u& u_param, Field_u& u_old_param):NumericalMethod(u_param,u_old_param){}
        virtual void solve(size_t j, double cfl) override{
           u[j] = 0.5*(((1+cfl)*u_old[j-1]) + ((1-cfl)*(u_old[j+1])));
        }
}; 
#endif