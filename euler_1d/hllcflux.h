#pragma once
#include <vector>
#include "eulerflux1d.h"
#include "conservedstate1d.h"
#include "eulerfluxcalculator.h"
#include "field.h"
#include "numericalflux.h"
class HLLCFlux : public NumericalFlux{
    protected:
        //double dx;
        double gamma;
    public:
        HLLCFlux(double gamma_):gamma(gamma_){}

        virtual EulerFlux1D numericalFlux(const ConservedState1D& UL, const ConservedState1D& UR,  double dt) override {
            
                EulerFlux1D fL = EulerFluxCalculator::fluxcalculate(UL, gamma);
                EulerFlux1D fR = EulerFluxCalculator::fluxcalculate(UR, gamma);
                //HLLC States Computation
                double ul = UL.velocity();
                double ur = UR.velocity();
                double pl = UL.pressure(gamma);
                double pr = UR.pressure(gamma);
                double rho_l = UL.rho;
                double rho_r = UR.rho;
                double al = UL.soundspeed(gamma);
                double ar = UR.soundspeed(gamma);
                double sl = std::min(ul-al, ur-ar);
                double sr = std::max(ul+al, ur+ar);

                double s_star = (pr - pl + (rho_l*ul*(sl-ul))- (rho_r*ur*(sr-ur)))/((rho_l*(sl-ul))- (rho_r*(sr-ur)) );
                double p_star = pl + (rho_l*(sl-ul)*(s_star-ul));

                double rho_star_l = rho_l *(sl-ul)/(sl-s_star);
                double rho_star_r = rho_r *(sr-ur)/(sr-s_star);

                ConservedState1D UL_star {
                    rho_star_l,
                    rho_star_l*s_star,
                    rho_star_l*((UL.rhoE/rho_l) + ((s_star-ul) * (s_star + (pl / (rho_l*(sl-ul))))) )
                };
                ConservedState1D UR_star {
                    rho_star_r,
                    rho_star_r*s_star,
                    rho_star_r*((UR.rhoE/rho_r) + ((s_star-ur) * (s_star + (pr / (rho_r*(sr-ur))))) )
                };

                EulerFlux1D FL_star = fL + (UL_star - UL) * sl;
                EulerFlux1D FR_star = fR + (UR_star - UR) * sr;
                EulerFlux1D F;

                if (sl >= 0){
                    F = fL;
                }
                else if (s_star >= 0 && sl <= 0){
                    F = FL_star;
                }
                else if (sr >= 0 && s_star <= 0){
                    F = FR_star;
                }else{
                    F = fR;
                }
                
                
                return F;
            
        }
};