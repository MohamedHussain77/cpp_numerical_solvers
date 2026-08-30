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
    
    virtual void step(Field& u,
                      std::vector<double>& rhs,
                      double dt) override { 
                        Field u_bar(u);
                        std::vector<double> rhs_2 (u.size());
                        for (size_t i{1}; i<(u.size()-1);++i){
                            u_bar[i] = u[i] + dt * rhs[i];} 
                        method.apply(u_bar, rhs_2, dt);
                        for (size_t i{1}; i<(u.size()-1);++i){
                            u[i] = u[i] + 0.5 * dt *(rhs[i]+rhs_2[i]);}                    
                      }
};
#endif