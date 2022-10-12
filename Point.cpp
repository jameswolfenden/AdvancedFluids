#include "Point.h"
#include <cmath>
#include "fluidConsts.h"

Point::Point(double p, double rho, double u)
{
    this->p = p;
    this->rho = rho;
    this->u = u;
    aCalc();
}

double Point::aCalc()
{
    a = sqrt((gamma * p) / rho);
    return a;
}
