#include "field.h"
#include "mesh.h"
#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H
class InitialCondition{
    protected:
        Field_u& u;
        Mesh& mesh;
    public:
        InitialCondition(Field_u& u_param, Mesh& mesh_param): u(u_param), mesh(mesh_param){}
        virtual void initialize() = 0;
};
#endif

