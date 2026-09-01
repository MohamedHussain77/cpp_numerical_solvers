#pragma once
#include "conservedstate1d.h"

struct EulerFlux1D{
    double massFlux;
    double momentumFlux;
    double energyFlux;

    EulerFlux1D operator+(const EulerFlux1D& other) const  {
        return{
        massFlux + other.massFlux,
        momentumFlux + other.momentumFlux,
        energyFlux + other.energyFlux};
    }
    EulerFlux1D operator-(const EulerFlux1D& other) const  {
        return{
        massFlux - other.massFlux,
        momentumFlux - other.momentumFlux,
        energyFlux - other.energyFlux};
    }
    EulerFlux1D operator-(const ConservedState1D& other) const {
        return{
        massFlux - other.rho,
        momentumFlux - other.rhou,
        energyFlux - other.rhoE};
    }
    EulerFlux1D operator*(double scalar) const {
        return{
        massFlux * scalar,
        momentumFlux * scalar,
        energyFlux * scalar};
    }
    EulerFlux1D operator/(double scalar) const {
        return{
        massFlux / scalar,
        momentumFlux / scalar,
        energyFlux / scalar};
    }
       
};
