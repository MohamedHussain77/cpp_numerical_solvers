#pragma once
#include <vector>
#include "eulerflux1d.h"
#include "conservedstate1d.h"
#include "eulerfluxcalculator.h"
#include "field.h"
#include "numericalflux.h"
class RoeFlux : public NumericalFlux{
    protected:
        //double dx;
        double gamma;
    public:
        RoeFlux(double gamma_):gamma(gamma_){}
        double HHcorrection (double lambda, double lambdaL, double lambdaR) const {
            double eps = std::max(0.0, std::max(lambda - lambdaL, lambdaR - lambda));
            if (std::abs(lambda) < eps)
                return (lambda*lambda + eps*eps) / (2.0*eps);
            return std::abs(lambda);
        }
        virtual EulerFlux1D numericalFlux(const ConservedState1D& UL, const ConservedState1D& UR,  double dt) override {
            
                EulerFlux1D fL = EulerFluxCalculator::fluxcalculate(UL, gamma);
                EulerFlux1D fR = EulerFluxCalculator::fluxcalculate(UR, gamma);
                //Roe States Computation
                double sl = std::sqrt(UL.rho);
                double sr = std::sqrt(UR.rho);
                double u_roe = (sl*UL.velocity()+sr*UR.velocity())/(sl+sr);
                double H_roe = (sl*UL.enthalpy(gamma) + sr*UR.enthalpy(gamma))/(sl+sr);
                double rho_roe = sl * sr;
                double a_roe = std::sqrt((gamma-1)*(H_roe - (0.5*u_roe*u_roe)));
                ConservedState1D r1_roe = {
                    1,
                    u_roe-a_roe,
                    H_roe-(u_roe*a_roe)
                };
                ConservedState1D r2_roe = {
                    1,
                    u_roe,
                    0.5*u_roe*u_roe
                };
                ConservedState1D r3_roe = {
                    1,
                    u_roe+a_roe,
                    H_roe+(u_roe*a_roe)
                };
                double rho_diff = UR.rho - UL.rho;
                double u_diff = UR.velocity() - UL.velocity();
                double p_diff = UR.pressure(gamma) - UL.pressure(gamma);
                double alpha1 = (p_diff - (rho_roe*a_roe*u_diff))/(2*a_roe*a_roe);
                double alpha2 = rho_diff - (p_diff/(a_roe*a_roe));
                double alpha3 = (p_diff + (rho_roe*a_roe*u_diff))/(2*a_roe*a_roe);
                /*ConservedState1D eigenvalue_L {
                    UL.velocity()-UL.soundspeed(gamma),
                    UL.velocity(),
                    UL.velocity()+UL.soundspeed(gamma)
                };
                ConservedState1D eigenvalue_R {
                    UR.velocity()-UR.soundspeed(gamma),
                    UR.velocity(),
                    UR.velocity()+UR.soundspeed(gamma)
                };
                ConservedState1D eigenvalue_roe {
                    u_roe-a_roe,
                    u_roe,
                    u_roe+a_roe
                };
                ConservedState1D delta = eigenvalue_R - eigenvalue_L;*/
                /*double delta1 = UR.velocity()-UR.soundspeed(gamma) - UL.velocity()+UL.soundspeed(gamma);
                double delta2 = UR.velocity() - UL.velocity();
                double delta3 = UR.velocity()+UR.soundspeed(gamma) - UL.velocity()-UL.soundspeed(gamma);*/
                

                double lambda1 = HHcorrection((u_roe-a_roe),UL.velocity()-UL.soundspeed(gamma),UR.velocity()-UR.soundspeed(gamma));
                double lambda2 = HHcorrection((u_roe),UL.velocity(),UR.velocity());
                double lambda3 = HHcorrection((u_roe+a_roe),UL.velocity()+UL.soundspeed(gamma),UR.velocity()+UR.soundspeed(gamma));


                

                ConservedState1D D =  (r1_roe * (alpha1 * lambda1)) 
                + (r2_roe * (alpha2 * lambda2))
                + (r3_roe * (alpha3 * lambda3)) ;


                EulerFlux1D F = (fL + fR) * 0.5 - (D * 0.5);
                return F;
            
        }
};

