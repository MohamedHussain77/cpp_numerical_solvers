#include "field.h"
#ifndef NUMERICALMETHOD_H
#define NUMERICALMETHOD_H

class NumericalMethod{
    protected:
        Field_u& u;
        Field_u& u_old;

    public:   
        NumericalMethod(Field_u& u_param, Field_u& u_old_param):u(u_param),u_old(u_old_param){}
        virtual void solve(size_t j, double cfl)=0;  
        


};
#endif
