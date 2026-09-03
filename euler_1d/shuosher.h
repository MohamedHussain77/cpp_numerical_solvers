#pragma once
#include "field.h"
#include "mesh.h"
#include <cmath>
#include "initialcondition.h"

class Shu_Osher_IC : public InitialCondition{
    //double rho;
    //double v;
    //double p;
    double gamma;
    public:
        Shu_Osher_IC(Field& u_, const Mesh& mesh_, double gamma_): InitialCondition(u_, mesh_),
        //rho(rho_), v(v_), p(p_),
        gamma(gamma_){}
        virtual void initialize() override {
            for(std::size_t i{0}; i < U.size(); ++i){
                double rho (1);
                double v (0);
                double p (1);
                if(mesh[i] < 0.1){
                    rho = 3.857143;
                    v = 2.629369;
                    p = 10.333333;
                }
                else{
                    rho = 1 + 0.2 * std::sin(50 * mesh[i] - 25) ;
                }
                U[i].rho = rho;
                U[i].rhou = rho * v;
                U[i].rhoE = (p/(gamma-1)) + (0.5* U[i].rhou * v);
            }
        }
};