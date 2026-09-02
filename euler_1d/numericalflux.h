#pragma once
#include <vector>
#include "eulerflux1d.h"
#include "conservedstate1d.h"
#include "eulerfluxcalculator.h"
#include "field.h"
class NumericalFlux{
    
    public:
        virtual EulerFlux1D numericalFlux(const ConservedState1D& UL, const ConservedState1D& UR, double dt) = 0;
        virtual ~NumericalFlux() = default;
};
