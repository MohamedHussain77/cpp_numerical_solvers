#include "field.h"
#ifndef NUMERICALMETHOD_H
#define NUMERICALMETHOD_H

class NumericalMethod{
    protected:
        Field_u& u;
        Field_u& u_old;
        double a;
        double delt;
        double delx;

    public:   
        NumericalMethod(Field_u& u_param, Field_u& u_old_param, double a_param, double delt_param, double delx_param):u(u_param),u_old(u_old_param), a(a_param),
        delt(delt_param), delx(delx_param){}
        virtual void solve()=0;  
        


};
#endif
