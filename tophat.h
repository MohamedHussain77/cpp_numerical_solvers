#include "initialcondition.h"
#ifndef TOPHAT_H
#define TOPHAT_H

class TopHat: public InitialCondition {
    private:
        double Amp; //Amplitude
        double xL; //left edge
        double xR; //right edge
    public:
        TopHat(Field_u& u_param, Mesh& mesh_param, double Amp_param, double xL_param, double xR_param):InitialCondition(u_param, mesh_param), Amp(Amp_param), 
        xL(xL_param), xR(xR_param){}

        virtual void initialize() override{
            for (size_t i{0}; i < u.size(); ++i){
                if (mesh[i]>=xL&&mesh[i]<=xR){
                    u[i]=Amp;
                }
                else {u[i]=0;}
            } 
        }
};
#endif