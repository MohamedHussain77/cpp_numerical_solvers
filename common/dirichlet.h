#include "field.h"
#include "boundarycondition.h"
#ifndef DIRICHLET_H
#define DIRICHLET_H
class Dirichlet: public BoundaryCondition{
    public:
        Dirichlet(Field& u_param, double C_param, BoundarySide side_param, double delx_param):BoundaryCondition(u_param,C_param,side_param, delx_param){}
        virtual void apply() override {
            if(side == BoundarySide::Left){
                u[0] = C;
                
            }else{
                u[u.size()-1] = C;
            }
                        
        }
};
#endif