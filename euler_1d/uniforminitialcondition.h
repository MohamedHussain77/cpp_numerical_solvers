#ifndef UNIFORMINITIALCONDITION_H
#define UNIFROMINITIALCONDITION_H
#include "field.h"
#include "mesh.h"
#include <cmath>
#include "initialcondition.h"

class UniformInitialCondition : public InitialCondition{
    double rho;
    double v;
    double p;
    double gamma;
    public:
        UniformInitialCondition(Field& u_, const Mesh& mesh_, double rho_, double v_, double p_, double gamma_): InitialCondition(u_, mesh_),
        rho(rho_), v(v_), p(p_), gamma(gamma_){}
        virtual void initialize() override {
            for(std::size_t i{0}; i < U.size(); ++i){
                U[i].rho = rho;
                U[i].rhou = rho * v;
                U[i].rhoE = (p/(gamma-1)) + (0.5* U[i].rhou * v);
            }
        }
};
#endif

