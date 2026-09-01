#pragma once
#include <cmath>
struct ConservedState1D{
    double rho;
    double rhou;
    double rhoE;
    double velocity() const {
            return rhou/rho;
        }
    double pressure(double gamma) const {
            return (gamma-1) * (rhoE - (0.5 * rhou * rhou/rho));
        }
    double soundspeed(double gamma) const {
        return std::sqrt((gamma*(gamma-1)/rho)*(rhoE - (0.5 * rhou * rhou / rho)));
    }
    ConservedState1D operator-(const ConservedState1D& other) const  {
        return{
        rho - other.rho,
        rhou - other.rhou,
        rhoE - other.rhoE};
    }
    ConservedState1D operator+(const ConservedState1D& other) const  {
        return{
        rho + other.rho,
        rhou + other.rhou,
        rhoE + other.rhoE};
    }
    ConservedState1D operator*(double scalar) const {
        return{
        rho * scalar,
        rhou * scalar,
        rhoE * scalar};
    }

};