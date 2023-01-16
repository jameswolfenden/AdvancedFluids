#include <cmath>
#include <iostream>
#include "2D/Domain2D.h"
#include "fluidConsts.h"
#include "2D/SolveRiemann.h"

int main()
{
    double rho, u, v, a, p, rhoL, uL, vL, aL, pL, rhoR, uR, vR, aR, pR;
    SolveRiemann rSolver;
    rhoL = 0.292967;
    rhoR = 2.7355e-11;
    pL = 12.8143;
    pR = 2.52242e-12;
    uL = -2.04943;
    uR = -44.8764;
    vL = 0;
    vR = 0;
    aL = sqrt((gammma * pL) / rhoL);
    aR = sqrt((gammma * pR) / rhoR);

    rSolver.findStar(rhoL, uL, vL, aL, pL, rhoR, uR, vR, aR,pR, rho, u, v, a, p);
    std::cout << "rho: " << rho << ", u: " << u << ", v: " << v << ", a: " << a << ", p: " << p << std::endl;
}