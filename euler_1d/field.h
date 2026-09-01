#pragma once
#include <vector>
#include <cstddef>
#include <utility>

#include "conservedstate1d.h"
#include "mesh.h"


class Field{
    private:
        std::vector<ConservedState1D> U;
        const Mesh& mesh;
        

    public:
        Field(const Mesh& mesh_param): mesh(mesh_param){
            U.resize(mesh.get_numcells());   
        }
        
        /*std::vector<double>& get_u(){
            return u;
        }
        const std::vector<double>& get_u() const {
            return u;
        }*/
        ConservedState1D& operator[](std::size_t j){
            return U[j];
        }
        const ConservedState1D& operator[](std::size_t j) const {
            return U[j];
        }
        void copy_from(const Field& u_param){
            U = u_param.U;
        }
        void swap_btw(Field& u_param){
            std::swap(U,u_param.U); //even though u_param.u is private member of different object, same class objects can access its private members
        }
        
        std::size_t size() const{
            return U.size();
        }
};

