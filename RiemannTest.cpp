#include <cmath>
#include <iostream>
#include "2D/Domain2D.h"
#include "fluidConsts.h"
#include "2D/SolveRiemann.h"

int main()
{
    double rho, u, v, a, p, dye, rhoL, uL, vL, aL, pL, rhoR, uR, vR, aR, pR;
    SolveRiemann rSolver;


    int caseTest = 33;

    if (caseTest == 1)
    {
        pL = 1.0;
        rhoL = 1.0;
        uL = 0.0;
        vL = 0.0;
        pR = 0.1;
        rhoR = 0.125;
        uR = 0.0;
        vR = 0.0;
    }
    else if (caseTest == 2)
    {
        pL = 0.4;
        rhoL = 1.0;
        uL = -2.0;
        vL = 0.0;
        pR = 0.4;
        rhoR = 1.0;
        uR = 2.0;
        vR = 0.0;
    }
    else if (caseTest == 3)
    {
        pL = 100.0;
        rhoL = 1.0;
        uL = 0.0;
        vL = 0.0;
        pR = 0.01;
        rhoR = 1.0;
        uR = 0.0;
        vR = 0.0;
    }
    else if (caseTest == 4)
    {
        pL = 0.01;
        rhoL = 1.0;
        uL = 0.0;
        vL = 0.0;
        pR = 100.0;
        rhoR = 1.0;
        uR = 0.0;
        vR = 0.0;
    }
    else if (caseTest == 5)
    {
        pL = 460.894;
        rhoL = 5.99924;
        uL = 19.5975;
        vL = 0.0;
        pR = 46.0950;
        rhoR = 5.99242;
        uR = -6.19633;
        vR = 0.0;
    }
    else {
        rhoL = 1.0;
        rhoR = 0.0;
        pL = 1.0;
        pR = 0.0;
        uL = 0.0;
        uR = 0.0;
        vL = 0.0;
        vR = 0.0;
    }

    if (rhoL==0.0)
        aL = 0.0;
    else
        aL = sqrt((gammma * pL) / rhoL);

    if (rhoR==0.0)
        aR = 0.0;
    else
        aR = sqrt((gammma * pL) / rhoL);

    rSolver.findStar(rhoL, uL, vL, aL, pL, 0.0, rhoR, uR, vR, aR, pR, 1.0, rho, u, v, a, p, dye);
    std::cout << "rho: " << rho << ", u: " << u << ", v: " << v << ", a: " << a << ", p: " << p << ", dye: " << dye << std::endl;
}