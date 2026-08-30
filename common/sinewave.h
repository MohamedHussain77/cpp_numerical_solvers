#include "initialcondition.h"
#include <cmath>
#ifndef SINEWAVE_H
#define SINEWAVE_H

class SineWave: public InitialCondition {
    private:
        // u(x,0) = Asin(kx+phi) + B
        double Amp; //Amplitude
        double k; //wavenumber 2pi/lambda
        double phi; //phase angle in radians
        double offset; //offset
    public:
        SineWave(Field& u_param, Mesh& mesh_param, double Amp_param, double k_param, double phi_param, double offset_param):InitialCondition(u_param, mesh_param), Amp(Amp_param), 
        k(k_param), phi(phi_param), offset(offset_param){}

        virtual void initialize() override{
            for (size_t i{0}; i < u.size(); ++i){
                u[i] = Amp * std::sin(k*mesh[i]+phi) + offset;
            } 
        }
};
#endif
