#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <numbers>
#include "mesh.h"
#include "field.h"
#include "solver.h"
#include "numericalmethod.h"
#include "upwind.h"
#include "laxwendroff.h"
#include "laxfriedrichs.h"
#include "dirichlet.h"
#include "neumann.h"
#include "robin.h"
#include "initialcondition.h"
#include "tophat.h"
#include "sinewave.h"
#include "gaussian.h"
#include "advection.h"
#include "diffusion.h"
#include "timestepper.h"
#include "expliciteulerstepper.h"
#include "diffusioncentraldiff.h"
int main(){
    double L = 1;
    int numcells = 200;
    Mesh mesh(numcells, L);

    double r = 0.4;
    double nu = 0.01;
    double totaltime = 0.5;
    
    Field_u u(mesh);
    GaussianPulse gauss(u, mesh, 1, 0.5, 0.05);
    gauss.initialize();
    //SineWave sine(u, mesh, 0.8, 4*std::numbers::pi,0,0);
    //sine.initialize();
    //Field_u u_analytical(mesh);
    //u_analytical.copy_from(u);
    //u.initialize_u_tophat();
    //Field_u u_old(mesh);
    Diffusioncentraldiff1D diffcentraldiff (u, mesh.get_delx(), nu);
    Dirichlet left(u, 0, BoundarySide::Left, mesh.get_delx());
    Dirichlet right(u, 0, BoundarySide::Right, mesh.get_delx()); //alpha * u + beta du/dx = C 
    ExplicitEulerStepper stepper;



    Diffusion solve(mesh, u, nu, r, totaltime, diffcentraldiff, stepper,  left, right);


    std::cout << "del x: " << mesh.get_delx() << std::endl; 
    //std::cout << "numsteps: " << solve.get_numsteps() << std::endl;

    solve.timeIntegrate();
    //sinewaveanalytical analyticsolution(u_analytical, mesh, totaltime, a, cfl, 0.8, 4*std::numbers::pi, 0,0);
    return 0; 
}
/*
class AnalyticalSolution_advection{
    protected:
        Field_u& u;
        Mesh& mesh;
        double totaltime;
        double a;
        double cfl;
        double delt;
        int numsteps;
    public:
        AnalyticalSolution_advection(Field_u& u_param, Mesh& mesh_param, double totaltime_param, double a_param, double cfl_param):
        u(u_param), mesh(mesh_param), totaltime(totaltime_param), a(a_param), cfl(cfl_param), delt(cfl*mesh.get_delx()/a), numsteps(static_cast<int>(totaltime/delt)){}
        virtual void analyticsln() = 0;
        
};

class sinewaveanalytical : public AnalyticalSolution_advection{
    private:
        // u(x,0) = Asin(kx+phi) + B
        double Amp; //Amplitude
        double k; //wavenumber 2pi/lambda
        double phi; //phase angle in radians
        double offset; //offset
    public: 
        sinewaveanalytical(Field_u& u_param, Mesh& mesh_param, double totaltime_param, double a_param, double cfl_param,
        double Amp_param, double k_param, double phi_param, double offset_param) : AnalyticalSolution_advection(u_param, mesh_param, totaltime_param, a_param, cfl_param), 
        Amp(Amp_param), k(k_param), phi(phi_param), offset(offset_param){};
        virtual void analyticsln() override {
            std::ofstream output("solution_ADVECTION_sine_analytical.csv");
            double time = 0.0;
            output << "t,x,u\n";
            for (size_t j{0}; j<mesh.get_numcells(); ++j){
                    output << time << "," << mesh[j] << "," << u[j] << "\n";
                }
            for (size_t i{0}; i<numsteps; ++i){
                //double time = delt * i;
                time += delt;
                for (size_t j{0}; j<mesh.get_numcells(); ++j){
                        u[j] = Amp * std::sin(k*(mesh[j]-(std::abs(a)*time))+phi) + offset;
                        output << time << "," << mesh[j] << "," << u[j] << "\n";

                    }  
            }
        }
};*/