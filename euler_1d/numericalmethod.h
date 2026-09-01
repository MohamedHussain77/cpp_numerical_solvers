#pragma once
#include "numericalflux.h"
#include "field.h"
#include "eulerflux1d.h"
#include "eulerfluxcalculator.h"
#include "boundarycondition.h"

class NumericalMethod{
    protected:
        //std::vector<EulerFlux1D>& F;
        //std::vector<ConservedState1D>& RHS;
        const Mesh& mesh;
        NumericalFlux& flux;
        BoundaryCondition& left;
        BoundaryCondition& right;
    
    public:
        NumericalMethod(const Mesh& mesh_, NumericalFlux& flux_, BoundaryCondition& left_, BoundaryCondition& right_):mesh(mesh_), flux(flux_), left(left_), right(right_){}
        //NumericalMethod(std::vector<EulerFlux1D>& F_, std::vector<EulerState1D>& RHS_):F(F_), RHS(RHS_){}
        void computeRHS(const Field& U, std::vector<EulerFlux1D>& F, std::vector<ConservedState1D>& RHS,  double dt){
            F.resize(U.size()+1);
            RHS.resize(U.size());
            ConservedState1D bcl = left.getBoundaryState(U);
            ConservedState1D bcr = right.getBoundaryState(U);
            for (std::size_t i{1}; i<U.size();++i){              
                F[i] = flux.numericalFlux(U[i-1], U[i], dt);
            }
            F[0] = flux.numericalFlux(bcl, U[0], dt);
            F[U.size()] = flux.numericalFlux(U[U.size()-1], bcr, dt);
            for (std::size_t i{0}; i<U.size();++i){
                EulerFlux1D fdiff = (F[i] - F[i+1])/mesh.get_delx();
                RHS[i].rho = fdiff.massFlux;
                RHS[i].rhou = fdiff.momentumFlux;
                RHS[i].rhoE = fdiff.energyFlux;


            }
        }
        
        
};