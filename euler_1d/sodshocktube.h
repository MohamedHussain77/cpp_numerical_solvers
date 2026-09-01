#pragma once
#include "field.h"
#include "mesh.h"
#include "initialcondition.h"
class SodShockTubeTest : public InitialCondition{
    double rho_l;
    double v_l;
    double p_l;
    double x;
    double rho_r;
    double v_r;
    double p_r;
    
    double gamma;
    public:
        SodShockTubeTest(Field& u_, const Mesh& mesh_, double rho_l_, double v_l_, double p_l_, double x_, double rho_r_, double v_r_, double p_r_, double gamma_): InitialCondition(u_, mesh_),
        rho_l(rho_l_), v_l(v_l_), p_l(p_l_), x(x_), rho_r(rho_r_), v_r(v_r_), p_r(p_r_), gamma(gamma_){}
        virtual void initialize() override {
            for(std::size_t i{0}; i < U.size(); ++i){
                double rho{0};
                double v{0};
                double p{1};
                if (mesh[i] <= x){
                    rho = rho_l;
                    v = v_l;
                    p = p_l;
                }else{
                    rho = rho_r;
                    v = v_r;
                    p = p_r;
                }
                U[i].rho = rho;
                U[i].rhou = rho * v;
                U[i].rhoE = (p/(gamma-1)) + (0.5* U[i].rhou * v);
            }
        }
};
