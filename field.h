#include <vector>
#include <cmath>
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
        std::size_t size() const{
            return u.size();
        }
};
#endif
