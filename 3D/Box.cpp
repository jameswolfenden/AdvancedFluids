#include "Box.h"
#include <cmath>
#include "../fluidConsts.h"
#include <iostream>

Box::Box(double p, double rho, double u, double v, double w)
{
    updatePrimatives(p, rho, u, v, w);
}
double Box::aCalc() // find speed of sound
{
    if (rho == 0.0)
        a = 0.0;
    else
        a = sqrt((gammma * p) / rho);
    return a;
}
// the following are to get the conservatives and fluxes from the primatives
double Box::u1()
{
    return rho;
}
double Box::u2()
{
    return rho * u;
}
double Box::u3()
{
    return rho * v;
}
double Box::u4()
{
    return rho * w;
}
double Box::u5()
{
    if (rho == 0.0)
        return 0.0;
    else
        return rho * (0.5 * (u * u + v * v + w * w) + p / ((gammma - 1) * rho));
}
double Box::f1()
{
    return rho * u;
}
double Box::f2()
{
    return rho * u * u + p;
}
double Box::f3()
{
    return rho * v * u;
}
double Box::f4()
{
    return rho * w * u;
}
double Box::f5()
{

    double result = u * (rho * (0.5 * (u * u + v * v * w * w) + p / ((gammma - 1) * rho)) + p);
    if (rho == 0.0)
        return 0.0;
    else

        if (result != result)
        std::cout << "u: " << u << ", rho: " << rho << ", v: " << v << ", p: " << p << std::endl;
    return result;
}
double Box::g1()
{
    return rho * v;
}
double Box::g2()
{
    return rho * v * u;
}
double Box::g3()
{
    return rho * v * v + p;
}
double Box::g4()
{
    return rho * w * v;
}
double Box::g5()
{
    if (rho == 0.0)
        return 0.0;
    else
        return v * (rho * (0.5 * (u * u + v * v + w * w) + p / ((gammma - 1) * rho)) + p);
}
double Box::h1()
{
    return rho * w;
}
double Box::h2()
{
    return rho * w * u;
}
double Box::h3()
{
    return rho * w * v;
}
double Box::h4()
{
    return rho * w * w + p;
}
double Box::h5()
{
    if (rho == 0.0)
        return 0.0;
    else
        return w * (rho * (0.5 * (u * u + v * v + w * w) + p / ((gammma - 1) * rho)) + p);
}

void Box::updateConservatives(double u1, double u2, double u3, double u4, double u5) // update the primatives from new conservative values
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
        double w_ = u4 / u1;
        double p_ = (gammma - 1) * (u5 - 0.5 * ((u2 * u2 + u3 * u3 + u4 * u4) / u1)); // these all need checking i guessed it ngl
        if (p_ < 0)
            std::cout << "p below zero, " << p_ << ", " << u_ << ", " << rho_ << ", " << u1 << ", " << u2 << ", " << u3 << ", " << u4 << std::endl;
        if (p_ != p_)
        {
            std::cout << "p error, " << p_ << ", " << u_ << ", " << rho_ << ", " << u1 << ", " << u2 << ", " << u3 << ", " << u4 << std::endl;
            updatePrimatives(p_, rho_, u_, v_, w_);
        }
        else
        {
            updatePrimatives(p_, rho_, u_, v_, w_);
        }
    }
}
void Box::updatePrimatives(double p, double rho, double u, double v, double w)
{
    this->u = u;
    this->v = v;
    this->w = w;
    this->p = p;
    this->rho = rho;
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
