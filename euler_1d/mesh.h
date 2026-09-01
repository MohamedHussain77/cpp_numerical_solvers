#include <vector>
#ifndef MESH_H
#define MESH_H
class Mesh{
    private:
        double L;
        int numcells;
        double delx ;
        std::vector<double> cellCenters;
        std::vector<double> faces;
        
    public:
        Mesh(int numcells_param, double L_param) : numcells(numcells_param), L(L_param), delx(L_param/numcells_param)
        {
            cellCenters.resize(numcells);
            faces.resize(numcells+1);
            for (std::size_t i{0}; i < numcells; ++i){
                cellCenters[i] = delx * (i+0.5);
                faces[i] = delx * i;
            }
            faces[numcells] = delx * numcells;
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
        const double& operator[](std::size_t j) const {
            return cellCenters[j];
        } 
        const double& operator()(std::size_t j) const {
            return faces[j];
        }    
};
#endif