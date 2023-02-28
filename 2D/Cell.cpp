#include "Cell.h"
#include <cmath>
#include "../fluidConsts.h"
#include <iostream>

Cell::Cell(double p, double rho, double u, double v, double dye)
{
    updatePrimatives(p, rho, u, v, dye);
}
double Cell::aCalc() // find speed of sound
{
    if (rho == 0.0)
        a = 0.0;
    else
        a = sqrt((gammma * p) / rho);
    return a;
}
// the following are to get the conservatives and fluxes from the primatives
double Cell::u1()
{
    return rho;
}
double Cell::u2()
{
    return rho * u;
}
double Cell::u3()
{
    return rho * v;
}
double Cell::u4()
{
    if (rho == 0.0)
        return 0.0;
    else
        return rho * (0.5 * (u * u + v * v) + p / ((gammma - 1) * rho));
}
double Cell::u5()
{
    return rho*dye;
}
double Cell::f1()
{
    return rho * u;
}
double Cell::f2()
{
    return rho * u * u + p;
}
double Cell::f3()
{
    return rho * v * u;
}
double Cell::f4()
{

    double result = u * (rho * (0.5 * (u * u + v * v) + p / ((gammma - 1) * rho)) + p);
    if (rho == 0.0)
        return 0.0;
    else

        if (result != result)
        std::cout << "u: " << u << ", rho: " << rho << ", v: " << v << ", p: " << p << std::endl;
    return result;
}
double Cell::f5()
{
    return rho*dye * u;
}
double Cell::g1()
{
    return rho * v;
}
double Cell::g2()
{
    return rho * v * u;
}
double Cell::g3()
{
    return rho * v * v + p;
}
double Cell::g4()
{
    if (rho == 0.0)
        return 0.0;
    else
        return v * (rho * (0.5 * (u * u + v * v) + p / ((gammma - 1) * rho)) + p);
}
double Cell::g5()
{
    return rho*dye * v;
}

void Cell::updateConservatives(double u1, double u2, double u3, double u4, double u5) // update the primatives from new conservative values
{
    if (u1 == 0)
    {
        updatePrimatives(0, 0, 0, 0, 0);
    }
    else
    {
        double rho_ = u1;
        double u_ = u2 / u1;
        double v_ = u3 / u1;
        double dye_ = u5 / u1;
        double p_ = (gammma - 1) * (u4 - 0.5 * ((u2 * u2 + u3 * u3) / u1)); // these all need checking i guessed it ngl
        if (p_ < 0)
            std::cout << "p below zero, " << p_ << ", " << u_ << ", " << rho_ << ", " << u1 << ", " << u2 << ", " << u3 << ", " << u4 << std::endl;
        if (p_ != p_)
        {
            std::cout << "p error, " << p_ << ", " << u_ << ", " << rho_ << ", " << u1 << ", " << u2 << ", " << u3 << ", " << u4 << std::endl;
            updatePrimatives(p_, rho_, u_, v_, dye_);
        }
        else
        {
            updatePrimatives(p_, rho_, u_, v_, dye_);
        }
    }
}
void Cell::updatePrimatives(double p, double rho, double u, double v, double dye)
{
    this->u = u;
    this->v = v;
    this->p = p;
    this->rho = rho;
    this->dye = dye;
    if (p >= 0.0 && rho >= 0.0)
    {
        aCalc();
        return;
    }
    if (p <= 0.0) // this is different to the number used in the vacuum shit to diagnose that the problem is not within the vacuum shit I CHANGED THIS
    {
        std::cout << "p is " << p << " rho is " << rho << std::endl;
        this->p = 0.0; // value not zero so other calcs still work I CHANGED THIS
        this->rho = 0.0;
    }
    if (rho <= 0.0)
    {
        std::cout << "p is " << p << " rho is " << rho << std::endl;
        this->rho = 0.0;
        this->p =0.0;
    }
    a=0.0;
    return;
}
