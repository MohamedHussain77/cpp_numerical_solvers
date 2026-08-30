#include "field.h"
#include "boundarycondition.h"
#ifndef NEUMANN_H
#define NEUMANN_H
class Neumann: public BoundaryCondition{
    public:
        Neumann(Field& u_param, double C_param, BoundarySide side_param, double delx_param):BoundaryCondition(u_param,C_param,side_param, delx_param){}
        virtual void apply() override {
            if(side == BoundarySide::Left){
                u[0] = u[1] - C * delx;
                
            }else{
                u[u.size()-1] = C * delx + u[u.size()-2] ;
            }
                        
        }
};
#endif


