#pragma once
#include "conservedstate1d.h"
class Limiter{
    public:
        Limiter() = default;
        virtual ConservedState1D sigma(const ConservedState1D& dU_minus, const ConservedState1D& dU_plus) = 0;
        double minmod(const double a, const double b)
        {
            if (a * b <= 0.0)
                return 0.0;

            return std::copysign(std::min(std::abs(a), std::abs(b)), a);
        }
        ConservedState1D minmod(const ConservedState1D& a,
                            const ConservedState1D& b)
            {
                return {
                    minmod(a.rho,  b.rho),
                    minmod(a.rhou, b.rhou),
                    minmod(a.rhoE, b.rhoE)
                };
            }
};
class VanLeer : public Limiter{
    public:
        VanLeer() = default;
        double vanLeer_scalar(const double dU_minus, const double dU_plus)
        {
            // Avoid division by zero
            if (dU_plus == 0.0)
                return 0.0;

            double r = dU_minus / dU_plus;

            double phi = (r + std::abs(r))
                    / (1.0 + std::abs(r));

            return phi * dU_plus;
        }
        virtual ConservedState1D sigma(const ConservedState1D& dU_minus, const ConservedState1D& dU_plus) override {
        return {
            vanLeer_scalar(dU_minus.rho,
                           dU_plus.rho),

            vanLeer_scalar(dU_minus.rhou,
                           dU_plus.rhou),

            vanLeer_scalar(dU_minus.rhoE,
                           dU_plus.rhoE)
        };
        }
};
class MinMod : public Limiter{
    public:
        MinMod() = default;
        virtual ConservedState1D sigma(const ConservedState1D& dU_minus, const ConservedState1D& dU_plus) override {
        ConservedState1D sigma_ = minmod(dU_minus, dU_plus); 
        return sigma_;
        }
};