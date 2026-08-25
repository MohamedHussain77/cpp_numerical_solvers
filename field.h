#include <vector>
#include "mesh.h"
#ifndef FIELD_H
#define FIELD_H
class Field_u{
    private:
        std::vector<double> u;
        const Mesh& mesh;

    public:
        Field_u(const Mesh& mesh_param): mesh(mesh_param){
            u.resize(mesh.get_numcells());   
        }
        void initialize_u_tophat(){
            for (size_t i{0}; i < mesh.get_numcells(); ++i){
                if (mesh[i]>=0.2&&mesh[i]<=0.4){
                    u[i]=1;
                }
                else {u[i]=0;}
            }     
        }
        /*std::vector<double>& get_u(){
            return u;
        }
        const std::vector<double>& get_u() const {
            return u;
        }*/
        double& operator[](size_t j){
            return u[j];
        }
        const double& operator[](size_t j) const {
            return u[j];
        }
        void copy_from(const Field_u& u_param){
            u = u_param.u;
        }
        void swap_btw(Field_u& u_param){
            std::swap(u,u_param.u); //even though u_param.u is private member of different object, same class objects can access its private members
        }
};
#endif