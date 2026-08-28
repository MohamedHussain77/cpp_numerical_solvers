#include "field.h"
#include "numericalmethod.h"
#ifndef LAXWENDROFF_H
#define LAXWENDROFF_H

class Laxwendroff: public NumericalMethod{

    public:
        Laxwendroff(Field_u& u_param, Field_u& u_old_param):NumericalMethod(u_param,u_old_param){}
        virtual void solve(size_t j, double cfl) override{
            u[j] = u_old[j] - ((0.5*cfl)*(u_old[j+1] - u_old[j-1]))+ ((0.5*cfl*cfl)*(u_old[j+1] - 2*u_old[j] + u_old[j-1]));
        }
}; 
#endif
