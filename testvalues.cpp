#include <cmath>
#include <iostream>
#include "fluidConsts.h"
#include "3D/SolveRiemann.h"

int main()
{
    double rho, u, v, w, a, p, rhoL, uL, vL, wL, aL, pL, rhoR, uR, vR, wR, aR, pR;
    SolveRiemann rSolver;

    for (int caseTest = 1; caseTest < 7; caseTest++)
    {
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
        else
        {
            rhoL = 1.0;
            rhoR = 0.0;
            pL = 1.0;
            pR = 0.0;
            uL = 0.0;
            uR = 0.0;
            vL = 0.0;
            vR = 0.0;
        }

        if (rhoL == 0.0)
            aL = 0.0;
        else
            aL = sqrt((gammma * pL) / rhoL);

        if (rhoR == 0.0)
            aR = 0.0;
        else
            aR = sqrt((gammma * pL) / rhoL);
        wL = 0.0;
        wR = 0.0;
        rSolver.findStar(rhoL, uL, vL, wL, aL, pL, rhoR, uR, vR, wR, aR, pR, rho, u, v, w, a, p);
        std::cout << "case: " << caseTest << ", rho: " << rho << ", u: " << u << ", v: " << v << ", a: " << a << ", p: " << p << std::endl;
    }
    return 0;
}