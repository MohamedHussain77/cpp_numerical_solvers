#pragma once
#include "conservedstate1d.h"
#include "limiter.h"
#include <utility>
#include <tuple>

class Reconstruction{
    protected:
        Limiter& limiter;
    public:
        Reconstruction(Limiter& limiter_): limiter(limiter_){}

        std::pair<ConservedState1D, ConservedState1D>
        construct(const ConservedState1D& ULL, const ConservedState1D& UL, const ConservedState1D& UR, const ConservedState1D& URR, double dx) {

            ConservedState1D ULnew = UL + limiter.sigma(UL-ULL, UR-UL) * 0.5;
            ConservedState1D URnew = UR - limiter.sigma(UR-UL, URR-UR) * 0.5;
            return {ULnew, URnew};

            
        }
        
};

