#include "Point.h"
#include <cmath>
#include "fluidConsts.h"

Point::Point(double p, double rho, double u)
{
    updatePrimatives(p,rho,u);
}
double Point::aCalc()
{
    a = sqrt((gamma * p) / rho);
    return a;
}
double Point::u1()
{
    return rho;
}
double Point::u2()
{
    return rho*u;
}
double Point::u3()
{
    return rho*(0.5*u*u+p/((gamma-1)*rho));
}
double Point::f1()
{
    return rho*u;
}
double Point::f2()
{
    return rho*u*u+p;
}
double Point::f3()
{
    return u*(rho*(0.5*u*u+p/((gamma-1)*rho))+p);
}
void Point::updateConservatives(double u1, double u2, double u3)
{
    rho = u1;
    u = u2/u1;
    p = (gamma-1)*(u3-0.5*((u2*u2)/u1));
    aCalc();
}
void Point::updatePrimatives(double p, double rho, double u)
{
    this->p = p;
    this->rho = rho;
    this->u = u;
    aCalc();
}