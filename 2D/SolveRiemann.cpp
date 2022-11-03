#include <cmath>
#include <vector>
#include "SolveRiemann.h"
#include "../fluidConsts.h"
#include <iostream>

double SolveRiemann::guessp(double &pL, double &rhoL, double &uL, double &vL, double &aL, double &pR, double &rhoR, double &uR, double &vR, double &aR)
{

    double QUSER = 2.0;

    double CUP = 0.25 * (rhoL + rhoR) * (aL + aR);
    double PPV = 0.5 * (pL + pR) + 0.5 * (uL - uR) * CUP;
    PPV = std::max(0.0, PPV);
    double PMIN = std::min(pL, pR);
    double PMAX = std::max(pL, pR);
    double QMAX = PMAX / PMIN;
    double PM = 0.0;

    if (QMAX < QUSER && (PMIN < PPV && PPV < PMAX))
    {
        // Select PVRS Riemann solver
        PM = PPV;
    }
    else if (PPV < PMIN)
    {
        // Select Two-Rarefaction Riemann solver
        double PQ = pow(pL / pR, G1);
        double UM = (PQ * uL / aL + uR / aR + G4 * (PQ - 1.0)) / (PQ / aL + 1 / aR);
        double PTL = 1 + G7 * (uL - UM) / aL;
        double PTR = 1 + G7 * (UM - uR) / aR;
        PM = 0.5 * (pL * pow(PTL, G3) + pR * pow(PTR, G3));
    }
    else
    {
        // Select Two-Shock Riemann solver with
        // PVRS as estimate
        double GEL = sqrt((G5 / rhoL) / (G6 * pL + PPV));
        double GER = sqrt((G5 / rhoR) / (G6 * pR + PPV));
        PM = (GEL * pL + GER * pR - (uR - uL)) / (GEL + GER);
    }

    return PM;
}

std::vector<double> SolveRiemann::prefun(double &P, double &PK, double &CK, double &DK)
{
    double F = 0.0;
    double FD = 0.0;
    if (P < PK)
    {
        double PRAT = P / PK;
        F = G4 * CK * (pow(PRAT, G1) - 1.0);
        FD = (1.0 / (DK * CK)) * pow(PRAT, (-G2));
    }
    else
    {

        double AK = G5 / DK;
        double BK = G6 * PK;
        double QRT = sqrt(AK / (BK + P));
        F = (P - PK) * QRT;
        FD = (1.0 - 0.5 * (P - PK) / (BK + P)) * QRT;
    }
    std::vector<double> fs = {F, FD};
    return fs;
}

double SolveRiemann::getRho(double &u, double &p, double &pL, double &rhoL, double &pR, double &rhoR)
{
    double rho;
    if (u >= 0) // pick rho value depending on the side of the discontinuity
    {
        if (p > pL)
            rho = rhoL * (((p / pL) + G6) / (G6 * (p / pL) + 1));
        else
            rho = rhoL * pow((p / pL), 1 / gammma);
    }
    else
    {
        if (p > pR)
            rho = rhoR * (((p / pR) + G6) / (G6 * (p / pR) + 1));
        else
            rho = rhoR * pow((p / pR), 1 / gammma);
    }
    return rho;
}

double getV(double u, double &uL, double &uR)
{
    double v;
    if (u >= 0) // pick rho value depending on the side of the discontinuity
        v = uL;
    else
        v = uR;
    return v;
}
std::vector<double> SolveRiemann::getStars(double &pL, double &rhoL, double &uL, double &vL, double &aL, double &pR, double &rhoR, double &uR, double &vR, double &aR)
{
    double POLD = guessp(pL, rhoL, uL, vL, aL, pR, rhoR, uR, vR, aR);
    std::cout << POLD << std::endl;
    double UDIFF = uR - uL;

    std::vector<double> FLs;
    std::vector<double> FRs;

    for (int i = 0; i < 100; i++) // change
    {
        FLs = prefun(POLD, pL, aL, rhoL);
        FRs = prefun(POLD, pR, aR, rhoR);

        double P = POLD - ((FLs[0] + FRs[0] + UDIFF) / (FLs[1] + FRs[1]));
        double CHANGE = 2.0 * std::abs((P - POLD) / (P + POLD));
  //      std::cout << "change " << CHANGE << std::endl;
  //      std::cout << "pold " << POLD << std::endl;
  //      std::cout << "p " << P << std::endl;
        POLD = P;
        if (CHANGE < (0.000001))
        {
  //          std::cout << "converge i=" << i << std::endl;
            break;
        }
    }
    double p = POLD;
    double u = 0.5 * (uL + uR) + 0.5 * (FRs[0] - FLs[0]);
    double rho = getRho(u, p, pL, rhoL, pR, rhoR);
    double v = getV(u, uL, uR);
    std::vector<double> results = {p, u, rho, v};
    return results;
}
