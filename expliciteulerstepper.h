#ifndef EXPLICITEULERSTEPPER_H
#define EXPLICITEULERSTEPPER_H
#include "timestepper.h"
#include "field.h"
#include <vector>
class ExplicitEulerStepper: public TimeStepper{
    public:
    virtual void step(Field_u& u,
                      std::vector<double>& rhs,
                      double dt) override {
                        for (size_t i{1}; i<(u.size()-1);++i){
                            u[i] += dt * rhs[i];} 

                      }

};
#endif