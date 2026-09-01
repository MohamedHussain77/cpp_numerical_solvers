#ifndef RK2STEPPER_H
#define RK2STEPPER_H
#include "timestepper.h"
#include "numericalmethod.h"
#include "field.h"
#include <vector>
class RK2Stepper: public TimeStepper{
    NumericalMethod& method;
    public:
        RK2Stepper(NumericalMethod& m) : method(m) {}
    
    virtual void step(Field& U,
                      std::vector<ConservedState1D>& rhs,
                      double dt) override { 
                        Field U_bar(U);
                        std::vector<ConservedState1D> rhs_2 (U.size());
                        std::vector<EulerFlux1D> F_bar (U.size()+1);
                        for (size_t i{1}; i<(U.size()-1);++i){
                            U_bar[i] = U[i] + (rhs[i]*dt);
                        } 
                        method.computeRHS(U_bar, F_bar, rhs_2, dt);
                        for (size_t i{1}; i<(U.size()-1);++i){
                            U[i] = U[i] + ((rhs[i]+rhs_2[i])*0.5 * dt);}                    
                      }
};
#endif