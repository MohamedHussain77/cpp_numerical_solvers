#include "field.h"
#include "boundarycondition.h"
#ifndef ROBIN_H
#define ROBIN_H
class Robin: public BoundaryCondition{
    protected:
        double alpha;
        double beta;
    public:
        Robin(Field& u_param, double C_param, BoundarySide side_param, double delx_param, double alpha_param, double beta_param):BoundaryCondition(u_param,C_param,side_param, delx_param), alpha(alpha_param), beta(beta_param){}
        virtual void apply() override{
            if(side == BoundarySide::Left){
                u[0] = (C - (beta * u[1]/delx))/(alpha - (beta/delx));
                
            }else{
                u[u.size()-1] = (C + (beta * u[u.size()-2]/delx))/(alpha + (beta/delx));
            }

        }
    
};
#endif