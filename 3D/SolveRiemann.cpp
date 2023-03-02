#include "SolveRiemann.h"
#include <cmath>
#include "../fluidConsts.h"
#include <iostream>

bool SolveRiemann::findStar(const double &rhoL, const double &uL, const double &vL, const double &wL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &wR, const double &aR, const double &pR, double &rho, double &u, double &v, double &w, double &a, double &p) // find the values at the faces between 2 cells adjacent in x
{
    if (testVacuum(rhoL, uL, vL, wL, aL, pL, rhoR, uR, vR, wR, aR, pR, rho, u, v, w, a, p))
        return true; // a vacuum is generated, values found without iteration required NEED TO CHECK SPEED OF SOUND!!!!!

    iterateP(rhoL, uL, vL, wL, aL, pL, rhoR, uR, vR, wR, aR, pR, rho, u, v, w, a, p);

    if (u >= 0.0) // pick rho value depending on the side of the discontinuity
    {
        if (p > pL)
            rho = rhoL * (((p / pL) + ((gammma - 1) / (gammma + 1))) / (((gammma - 1) / (gammma + 1)) * (p / pL) + 1));
        else
            rho = rhoL * pow((p / pL), 1 / gammma);
        v = vL;
        w = wL;
    }
    else
    {
        if (p > pR)
            rho = rhoR * (((p / pR) + ((gammma - 1) / (gammma + 1))) / (((gammma - 1) / (gammma + 1)) * (p / pR) + 1));
        else
            rho = rhoR * pow((p / pR), 1 / gammma);
        v = vR;
        w = wR;
    }
    a = sqrt((gammma * p) / rho);
    if (rho != rho)
        std::cout << "rhoerror u: " << u << ", rhoL: " << rhoL << ", rhoR: " << rhoR << ", vL: " << vL << ", vR: " << vR << ", pL: " << pL << ", pR: " << pR << ", p:" << p << std::endl;
    return false;
}

bool SolveRiemann::testVacuum(const double &rhoL, const double &uL, const double &vL, const double &wL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &wR, const double &aR, const double &pR, double &rho, double &u, double &v, double &w, double &a, double &p)
{
    if ((rhoL == 0.0 || pL == 0.0) || (rhoR == 0.0 || pR == 0.0)) //|| (2 * aL / (gammma - 1) + 2 * aR / (gammma - 1)) <= (uR - uL))
    {
        if (((rhoL == 0.0 || pL == 0.0) && (rhoR != 0.0 || pR != 0.0)) || (((rhoR != 0.0 || pR != 0.0) && (rhoL != 0.0 || pL != 0.0)) && 0 >= (uR - 2 * aR / (gammma - 1)))) // vacuum left not right
        {
            if (0 >= (uR + aR))
            {
                std::cout << "W_R" << std::endl;
                rho = rhoR;
                u = uR;
                v = vR; // made up
                w = wR; // made up
                p = pR;
                a = sqrt((gammma * p) / rho);
            }
            else if (0 <= (uR - 2 * aR / (gammma - 1)))
            {
                // std::cout << "W_0" << std::endl;
                // rho = 0.0;
                // p = 0.0;
                // u = 0.5 * (uL + uR); // SOURCE: I MADE IT UP NEED TO ASK
                // v = 0.5 * (vL + vR); // SOURCE: I MADE IT UP NEED TO ASK
                // u = uL;              // idk
                // a = 0.0;
                std::cout << "W_L" << std::endl;
                rho = rhoL;
                p = pL;
                u = uL;
                v = vL;
                w = wL;
                a = aL;
            }
            else
            {
                std::cout << "W_RFan" << std::endl;
                rho = rhoR * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (uR) / aR, 2 / (gammma - 1));
                u = 2 / (gammma + 1) * (-aR + ((gammma - 1) / 2) * uR);
                v = vR * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (uR) / aR, 2 / (gammma - 1));
                w = wR * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (uR) / aR, 2 / (gammma - 1));
                p = pR * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (uR) / aR, (2 * gammma) / (gammma - 1));
                a = sqrt((gammma * p) / rho);
            }
            return true;
        }
        else if (((rhoR == 0.0 || pR == 0.0) && (rhoL != 0.0 || pL != 0.0)) || (((rhoR != 0.0 || pR != 0.0) && (rhoL != 0.0 || pL != 0.0)) && 0 <= (uL + 2 * aL / (gammma - 1)))) // vacuum right not left
        {
            if (0 <= (uL - aL))
            {
                std::cout << "W_L" << std::endl;
                rho = rhoL;
                u = uL;
                v = vL; // made up
                w = wL; // made up
                p = pL;
                a = sqrt((gammma * p) / rho);
            }
            else if (0 >= (uL + 2 * aL / (gammma - 1)))
            {
                // std::cout << "W_0" << std::endl;
                // rho = 0.0;
                // p = 0.0;
                // u = 0.5 * (uL + uR); // SOURCE: I MADE IT UP NEED TO ASK
                // v = 0.5 * (vL + vR); // SOURCE: I MADE IT UP NEED TO ASK
                // u = uR;              // idk
                // a = 0.0;
                std::cout << "W_R" << std::endl;
                rho = rhoR;
                p = pR;
                u = uR;
                v = vR;
                w = wR;
                a = aR;
            }
            else
            {
                std::cout << "W_LFan" << std::endl;
                rho = rhoL * pow(2 / (gammma + 1) + (gammma - 1) / (gammma + 1) * (uL) / aL, 2 / (gammma - 1));
                u = 2 / (gammma + 1) * (aL + ((gammma - 1) / 2) * uL);
                v = vL * pow(2 / (gammma + 1) + (gammma - 1) / (gammma + 1) * (uL) / aL, 2 / (gammma - 1));
                w = wL * pow(2 / (gammma + 1) + (gammma - 1) / (gammma + 1) * (uL) / aL, 2 / (gammma - 1));
                p = pL * pow(2 / (gammma + 1) + (gammma - 1) / (gammma + 1) * (uL) / aL, (2 * gammma) / (gammma - 1));
                a = sqrt((gammma * p) / rho);
            }
            return true;
        }
        else if (((rhoR == 0.0 || pR == 0.0) && (rhoL == 0.0 || pL == 0.0)) || ((uL + 2 * aL / (gammma - 1)) < 0 && (uR - 2 * aR / (gammma - 1)) > 0)) // vacuum left and right
        {
            std::cout << "W_0" << std::endl;
            rho = 0.0;
            p = 0.0;
            u = 0.5 * (uL + uR); // SOURCE: I MADE IT UP NEED TO ASK
            v = 0.5 * (vL + vR); // SOURCE: I MADE IT UP NEED TO ASK
            w = 0.5 * (wL + wR); // SOURCE: I MADE IT UP NEED TO ASK
            a = 0.0;
            return true;
        }
    }
    return false;
}

bool SolveRiemann::pickStartVal(const int errorStage, const double &rhoL, const double &uL, const double &vL, const double &wL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &wR, const double &aR, const double &pR, double &rho, double &u, double &v, double &w, double &a, double &p)
{

    double G1 = (gammma - 1) / (2 * gammma);
    double G2 = (gammma + 1) / (2 * gammma);
    double G3 = 2 * gammma / (gammma - 1);
    double G4 = 2 / (gammma - 1);
    double G5 = 2 / (gammma + 1);
    double G6 = (gammma - 1) / (gammma + 1);
    double G7 = (gammma - 1) / 2;
    double G8 = gammma - 1;

    double q_user = 2.0;

    double c_up = 0.25 * (rhoL + rhoR) * (aL + aR);
    double p_PV = 0.5 * (pL + pR) + 0.5 * (uL - uR) * c_up;
    p_PV = std::max(0.0, p_PV);
    double p_min = std::min(pL, pR);
    double p_max = std::max(pL, pR);
    double q_max = p_max / p_min;

    // different guesses for p
    if (errorStage == 0)
    {

        if (q_max < q_user && (p_min < p_PV && p_PV < p_max))
        {
            // Select PVRS Riemann solver
            p = p_PV;
        }
        else if (p_PV < p_min)
        {
            // Select Two-Rarefaction Riemann solver
            double p_q = pow(pL / pR, G1);
            double u_m = (p_q * uL / aL + uR / aR + G4 * (p_q - 1.0)) / (p_q / aL + 1 / aR);
            double p_TL = 1 + G7 * (uL - u_m) / aL;
            double p_TR = 1 + G7 * (u_m - uR) / aR;
            p = 0.5 * (pL * pow(p_TL, G3) + pR * pow(p_TR, G3));
        }
        else
        {
            // Select Two-Shock Riemann solver with
            // PVRS as estimate
            double ge_L = sqrt((G5 / rhoL) / (G6 * pL + p_PV));
            double ge_R = sqrt((G5 / rhoR) / (G6 * pR + p_PV));
            p = (ge_L * pL + ge_R * pR - (uR - uL)) / (ge_L + ge_R);
        }
    }
    else if (errorStage == 1)
    {
        std::cout << "went to error stage 1" << std::endl;
        p = pow((aL + aR - 0.5 * (gammma - 1) * (uR - uL)) / (aL / pow(pL, (gammma - 1) / (2 * gammma)) + aR / pow(pR, (gammma - 1) / (2 * gammma))), (2 * gammma) / (gammma - 1));
    }
    else if (errorStage == 2)
    {
        p = 0.5 * (pL + pR);
    }
    else if (errorStage == 3)
    {
        double p_PV = 0.5 * (pL + pR) + 0.5 * (uL - uR) * 0.5 * (rhoL + rhoR) * 0.5 * (aL + aR);
        // if (p_PV > TOL)
        p = p_PV;
        //  else
        //     p = TOL;
    }
    else if (errorStage == 4)
    {
        double p_PV = 0.5 * (pL + pR) + 0.5 * (uL - uR) * 0.5 * (rhoL + rhoR) * 0.5 * (aL + aR);
        if (p_PV > TOL)
            p = p_PV;
        else
            p = TOL;
        double AL = 2 / ((gammma + 1) * rhoL);
        double BL = pL * (gammma - 1) / (gammma + 1);
        double AR = 2 / ((gammma + 1) * rhoR);
        double BR = pR * (gammma - 1) / (gammma + 1);
        double gL = pow(AL / (p + BL), 0.5);
        double gR = pow(AR / (p + BR), 0.5);
        double p_TS = (gL * pL + gR * pR - (uR - uL)) / (gL + gR);
        //   if (p_TS > TOL)
        p = p_TS;
        //  else
        //      p = TOL;
    }
    else if (errorStage == 5)
    {
        p = p_PV;
    }
    else if (errorStage == 6)
    {
        // Select Two-Rarefaction Riemann solver
        double p_q = pow(pL / pR, G1);
        double u_m = (p_q * uL / aL + uR / aR + G4 * (p_q - 1.0)) / (p_q / aL + 1 / aR);
        double p_TL = 1 + G7 * (uL - u_m) / aL;
        double p_TR = 1 + G7 * (u_m - uR) / aR;
        p = 0.5 * (pL * pow(p_TL, G3) + pR * pow(p_TR, G3));
    }
    else if (errorStage == 7)
    {
        // Select Two-Shock Riemann solver with
        // PVRS as estimate
        double ge_L = sqrt((G5 / rhoL) / (G6 * pL + p_PV));
        double ge_R = sqrt((G5 / rhoR) / (G6 * pR + p_PV));
        p = (ge_L * pL + ge_R * pR - (uR - uL)) / (ge_L + ge_R);
    }
    else if (errorStage == 8)
    {
        p = 1 / (rhoL * aL + rhoR * aR) * (rhoR * aR * pL + rhoL * aL * pR + rhoL * aL * rhoR * aR * (uL - aR));
    }
    else if (errorStage == 9)
    {
        p = 1 * pow(10, -6);
    }
    else
    {
        std::cout << "Error converging on p in x" << std::endl;
        std::cout << "pL: " << pL << ", rhoL: " << rhoL << ", uL: " << uL << ", aL: " << aL << ", pR: " << pR << ", rhoR: " << rhoR << ", uR: " << uR << ", aR: " << aR << std::endl;
        return false; // this doesnt actually stop the program but if you see that in the console the timestep is probably too small
    }
    return true;
}
void SolveRiemann::iterateP(const double &rhoL, const double &uL, const double &vL, const double &wL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &wR, const double &aR, const double &pR, double &rho, double &u, double &v, double &w, double &a, double &p)
{
    bool iterate = true;
    double fsL, fsR, d_fsL, d_fsR, change;
    int errorStage = 0;
    int count = 0;

    while (iterate) // loop to try all the initial p values
    {
        pickStartVal(errorStage, rhoL, uL, vL, wL, aL, pL, rhoR, uR, vR, wR, aR, pR, rho, u, v, w, a, p);
        count = 0;
        // std::cout << "p guessed: " << p<< std::endl;
        while (iterate) // loop to iterate the p value
        {
            errno = 0;
            double A = 2 / ((gammma + 1) * rhoL);
            double B = pL * (gammma - 1) / (gammma + 1);
            if (p > pL) // shock
            {
                fsL = (p - pL) * pow(A / (p + B), 0.5);
                d_fsL = pow(A / (p + B), 0.5) * (1 - (p - pL) / (2 * (B + p)));
            }
            else // expansion
            {
                fsL = 2 * aL / (gammma - 1) * (pow(p / pL, (gammma - 1) / (2 * gammma)) - 1);
                d_fsL = 1 / (rhoL * aL) * pow(p / pL, -(gammma + 1) / (2 * gammma));
            }
            A = 2 / ((gammma + 1) * rhoR);
            B = pR * (gammma - 1) / (gammma + 1);
            if (p > pR) // shock
            {
                fsR = (p - pR) * pow(A / (p + B), 0.5);
                d_fsR = pow(A / (p + B), 0.5) * (1 - (p - pR) / (2 * (B + p)));
            }
            else // expansion
            {
                fsR = 2 * aR / (gammma - 1) * (pow(p / pR, (gammma - 1) / (2 * gammma)) - 1);
                d_fsR = 1 / (rhoR * aR) * pow(p / pR, -(gammma + 1) / (2 * gammma));
            }
            if (errno != 0) // p is probably negative giving an issue with pow, we need the next guess of p
            {
                errorStage++;
                break; // exit the iterate loop to try next p
            }

            double f = fsL + fsR - uL + uR;
            double d_f = d_fsL + d_fsR;
            change = f / d_f;
            //    std::cout << f << ", " << d_f << ", " << change << std::endl;
            p = p - change; // Update new estimate of p*
            count++;

            if (TOL >= 2 * fabs(change / (change + 2 * p))) // iteration limit (slightly different to notes as abs of entire rhs)
            {
                iterate = false; // we have converged
            }
            if (count > 10000)
            {
                errorStage++;
                break;
            }
        }
    }
    u = 0.5 * (uL + uR) + 0.5 * (fsR - fsL); // u*
    return;
}