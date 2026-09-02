#pragma once
#include <vector>
#include "eulerflux1d.h"
#include "conservedstate1d.h"
#include "eulerfluxcalculator.h"
#include "field.h"
#include "numericalflux.h"
class LaxFriedrichsFlux : public NumericalFlux{
    protected:
        double dx;
        double gamma;
    public:
        LaxFriedrichsFlux(double dx_, double gamma_):dx(dx_),gamma(gamma_){}
        virtual EulerFlux1D numericalFlux(const ConservedState1D& UL, const ConservedState1D& UR,  double dt) override {
            //for(std::size_t i{1}; i<F.size()-1; ++i){
                EulerFlux1D fL = EulerFluxCalculator::fluxcalculate(UL, gamma);
                EulerFlux1D fR = EulerFluxCalculator::fluxcalculate(UR, gamma);
                EulerFlux1D F = (fL + fR) * 0.5 - ((UR - UL) * (0.5*dx/dt));
                return F;
            //}
            //F[0] = ((EulerFluxCalculator::fluxcalculate(ghostL, gamma) 
            //+ EulerFluxCalculator::fluxcalculate(U[0], gamma)) * 0.5) 
            //- ((U[0] - ghostL) * (0.5*dx/dt));
            //F[F.size()-1] = ((EulerFluxCalculator::fluxcalculate(ghostR, gamma) 
            //+ EulerFluxCalculator::fluxcalculate(U[U.size()-1], gamma)) * 0.5) 
            //- ((ghostR - U[U.size()-1]) * (0.5*dx/dt));
        }
};
