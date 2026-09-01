#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H
#include "field.h"
#include "mesh.h"

class InitialCondition{
    protected:
        Field& U;
        const Mesh& mesh;
    public:
        InitialCondition(Field& u_, const Mesh& mesh_): U(u_), mesh(mesh_){}
        virtual void initialize() = 0;
};
#endif


