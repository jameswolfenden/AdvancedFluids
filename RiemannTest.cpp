#include <cmath>
#include <iostream>
#include "2D/Domain2D.h"
#include "fluidConsts.h"
#include "2D/SolveRiemann.h"

int main()
{
    double rho, u, v, a, p, rhoL, uL, vL, aL, pL, rhoR, uR, vR, aR, pR;
    SolveRiemann rSolver;
    rhoL = 0.00843769;
    rhoR = 1.0772e-08;
    pL = 4.63201e-07;
    pR = 7.15694e-12;
    uL = 1.01693;
    uR = -0.0564724;
    vL = 0;
    vR = 0;
    aL = sqrt((gammma * pL) / rhoL);
    aR = sqrt((gammma * pR) / rhoR);

    rSolver.findStar(rhoL, uL, vL, aL, pL, rhoR, uR, vR, aR,pR, rho, u, v, a, p);
    std::cout << "rho: " << rho << ", u: " << u << ", v: " << v << ", a: " << a << ", p: " << p << std::endl;
}