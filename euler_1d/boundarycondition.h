#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H
#include "field.h"
#include <cmath>

enum class BoundarySide {
    Left,
    Right
};

class BoundaryCondition{
    protected:
        //Field_Euler1D& u;
        //EulerState1D& bc;
        //double C1;
        //double C2;
        //double C3;
        BoundarySide side;
        //ConservedState1D ghostcell;
        //double dx;
        
        

    public:
        BoundaryCondition(BoundarySide side_param): side(side_param){}
        virtual ConservedState1D getBoundaryState(const Field& U) const =0;

};


class zeroGradient : public BoundaryCondition{
    public:
        zeroGradient(BoundarySide side_param): BoundaryCondition(side_param){}
        virtual ConservedState1D getBoundaryState(const Field& U) const override {
            if (side==BoundarySide::Left){
                return U[0];

            }else{
                return U[U.size()-1];
            }      

        }

};
class periodic : public BoundaryCondition{
    public:
        periodic(BoundarySide side_param): BoundaryCondition(side_param){}
        virtual ConservedState1D getBoundaryState(const Field& U) const override {
            if (side==BoundarySide::Left){
                return U[U.size()-1];

            }else{
                return U[0];
            }      

        }

};
class slipWall : public BoundaryCondition{
    public:
        slipWall(BoundarySide side_param): BoundaryCondition(side_param){}
        virtual ConservedState1D getBoundaryState(const Field& U) const override {
            if (side==BoundarySide::Left){
                return {
                    U[0].rho, -U[0].rhou, U[0].rhoE
                };

            }else{
                return {
                    U[U.size()-1].rho, -U[U.size()-1].rhou, U[U.size()-1].rhoE
                };
            }      

        }

};
class SupersonicInletBC : public BoundaryCondition{
    protected:
        double rho;
        double u;
        double p;
        double gamma;

    public:
        SupersonicInletBC(BoundarySide side_param, double rho_, double u_, double p_, double gamma_): BoundaryCondition(side_param),
        rho(rho_), u(u_), p(p_), gamma(gamma_){}
        virtual ConservedState1D getBoundaryState(const Field& U) const override {
            if (side==BoundarySide::Left){
                return {
                    rho, rho*u, p/(gamma-1)+(0.5*rho*u*u)
                };

            }else{
                return {
                    rho, -rho*u, p/(gamma-1)+(0.5*rho*u*u)
                };
            }      

        }

};
class SubsonicInletBC : public BoundaryCondition{
    protected:
        //double rho;
        double u;
        double p;
        double gamma;

    public:
        SubsonicInletBC(BoundarySide side_param, double u_, double p_, double gamma_): BoundaryCondition(side_param),
        u(u_), p(p_), gamma(gamma_){}
        virtual ConservedState1D getBoundaryState(const Field& U) const override {
            if (side==BoundarySide::Left){
                //double u0 = ;
                //double p_out = U[0].pressure(gamma);
                double c0 = std::sqrt(gamma*U[0].pressure(gamma)/U[0].rho);
                double cin = 0.5 * (gamma - 1) * (u - U[0].velocity() + (2*c0/(gamma-1)));
                double rhoin = gamma*p/(cin*cin);
                return {
                    rhoin , rhoin *u, p/(gamma-1)+(0.5*rhoin*u*u)
                };

            }else{
                double c0 = std::sqrt(gamma*U[U.size()-1].pressure(gamma)/U[U.size()-1].rho);
                double cin = 0.5 * (gamma - 1) * (U[U.size()-1].velocity()+(2*c0/(gamma-1))+u);
                double rhoin = gamma*p/(cin*cin);
                return {
                    rhoin , -rhoin *u, p/(gamma-1)+(0.5*rhoin*u*u)
                };
            }      

        }

};

class SubsonicOutletBC : public BoundaryCondition{
    protected:
        double p_out;
        double gamma;
    public:
        SubsonicOutletBC(BoundarySide side_param, double p_out_, double gamma_): BoundaryCondition(side_param), p_out(p_out_), gamma(gamma_){}
        virtual ConservedState1D getBoundaryState(const Field& U) const override {
            if (side==BoundarySide::Left){
                double K = U[0].pressure(gamma)/ std::pow(U[0].rho, gamma);
                double rhoout = std::pow(p_out/K, (1/gamma));
                double cout = std::sqrt(gamma*p_out/rhoout);
                double uout = +(2*cout/(gamma-1)) + U[0].velocity() - ((2/(gamma-1)) * std::sqrt(gamma*U[0].pressure(gamma)/U[0].rho)) ;
                return {
                    rhoout, rhoout*uout, p_out/(gamma-1)+(0.5*rhoout*uout*uout)                  
                };

            }else{
                double K = U[U.size()-1].pressure(gamma)/ std::pow(U[U.size()-1].rho, gamma);
                double rhoout = std::pow(p_out/K, (1/gamma));
                double cout = std::sqrt(gamma*p_out/rhoout);
                double uout = -(2*cout/(gamma-1)) + U[U.size()-1].velocity() + ((2/(gamma-1)) * std::sqrt(gamma*U[U.size()-1].pressure(gamma)/U[U.size()-1].rho)) ;
                return {
                    rhoout, rhoout*uout, p_out/(gamma-1)+(0.5*rhoout*uout*uout)
                };
            }      

        }

};
#endif


