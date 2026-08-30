#include <vector>
#ifndef MESH_H
#define MESH_H
class Mesh{
    private:
        double L;
        int numcells;
        double delx ;
        std::vector<double> x;
        
    public:
        Mesh(int numcells_param, double L_param) : numcells(numcells_param), L(L_param), delx(L_param/numcells_param)
        {
            x.resize(numcells);
            for (size_t i{0}; i < numcells; ++i){
                x[i] = delx * (i+0.5);
            }
        }
        int get_numcells() const{
            return numcells;    
        }
        double get_L() const{
            return L;
        }
        double get_delx() const{
            return delx;
        }
        /*const std::vector<double>& get_x() const{
            return x;
        }*/
        const double& operator[](size_t j) const {
            return x[j];
        }     
};
#endif