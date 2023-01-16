#include <cmath>
#include <iostream>
#include "Domain2D.h"
#include "../fluidConsts.h"
#include "SolveRiemann.h"

int main()
{
    double rho, u, v, a, p;
    SolveRiemann rSolver;
    rSolver.findStar(1.0, 0.0, 0.0, sqrt((gammma * 1.0) / 1.0), 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, rho, u, v, a, p);
    std::cout << "rho: " << rho << ", u: " << u << ", v: " << v << ", a: " << a << ", p: " << p << std::endl;
}