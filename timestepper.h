#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H
#include "field.h"
#include <vector>
class TimeStepper{
    public:
    virtual ~TimeStepper() = default;
    virtual void step(Field_u& u,
                      std::vector<double>& rhs,
                      double dt) = 0;
};
#endif
