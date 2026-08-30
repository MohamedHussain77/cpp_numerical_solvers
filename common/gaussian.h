#include "initialcondition.h"
#include <cmath>
#ifndef GAUSSIAN_H
#define GAUSSIAN_H

class GaussianPulse: public InitialCondition {
    private:
        // u(x,0) = A exp [- (x-xc)^2/(2sigma^2)]
        double Amp; //Amplitude
        double xc; //center of gaussian
        double w; //width
        
    public:
        GaussianPulse(Field& u_param, Mesh& mesh_param, double Amp_param, double xc_param, double w_param):InitialCondition(u_param, mesh_param), Amp(Amp_param), 
        xc(xc_param), w(w_param){}

        virtual void initialize() override{
            for (size_t i{0}; i < u.size(); ++i){
                u[i] = Amp * std::exp(-0.5*std::pow(((mesh[i]-xc)/w),2));
            } 
        }
};
#endif