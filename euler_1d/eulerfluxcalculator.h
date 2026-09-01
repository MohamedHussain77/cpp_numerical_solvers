#pragma once
#include "eulerflux1d.h"
#include "field.h"

class EulerFluxCalculator {
    public:
        
        static EulerFlux1D fluxcalculate(const ConservedState1D& U, double gamma ){
                return{ U.rhou,
                U.rho*U.velocity()*U.velocity()+U.pressure(gamma),
                U.velocity()*(U.rhoE+U.pressure(gamma))};
        };
};