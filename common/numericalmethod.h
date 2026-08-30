#ifndef NUMERICALMETHOD_H
#define NUMERICALMETHOD_H

#include "field.h"
#include <vector>

class NumericalMethod {
protected:
    //Field& u;
    //std::vector<double>& rhs;
    double dx;

public:
    NumericalMethod(//Field& u_param,
                    //std::vector<double>& rhs_param,
                    double dx_param)
        :  dx(dx_param) {}

    virtual ~NumericalMethod() = default;

    // Compute RHS = L(u) at current time level.
    // dt may be needed for some schemes (e.g. Lax–Wendroff).
    virtual void apply(Field& u, std::vector<double>& rhs, double dt) = 0;
};

#endif