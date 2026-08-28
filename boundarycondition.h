#include "field.h"
#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H
enum class BoundarySide {
    Left,
    Right
};

class BoundaryCondition{
    protected:
        Field_u& u;
        double C;
        BoundarySide side;
        double delx;
        

    public:
        BoundaryCondition(Field_u& u_param, double C_param, BoundarySide side_param, double delx_param): u(u_param), C(C_param), side(side_param), delx(delx_param){}
        virtual void apply()=0;

};
#endif