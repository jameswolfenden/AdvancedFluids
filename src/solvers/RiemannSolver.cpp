#include "solvers/RiemannSolver.hpp"
#include <cmath>
#include <iostream>

namespace fluid
{

    bool RiemannSolver::findStar(const double &rhoL, const double &uL, const double &vL, const double &wL,
                                 const double &aL, const double &pL, const double &rhoR, const double &uR,
                                 const double &vR, const double &wR, const double &aR, const double &pR,
                                 Flux &fl)
    {
        if (testVacuum(rhoL, uL, vL, wL, aL, pL, rhoR, uR, vR, wR, aR, pR, fl))
        {
            std::cout << "vacuum" << std::endl;
            return true; // a vacuum is generated, values found without iteration required
        }

        double tempP;
        double tempU;
        if (!iterateP(rhoL, uL, aL, pL, rhoR, uR, aR, pR, tempP, tempU))
        {
            std::cout << "iteration failed" << std::endl;
            return false; // iteration failed, abort
        }

        if (!pickSide(rhoL, vL, wL, pL, rhoR, vR, wR, pR, fl, tempP, tempU))
        {
            std::cout << "pick side failed, rho is nan" << std::endl;
            return false; // rho is nan
        }
        return true;
    }

    bool RiemannSolver::pickSide(const double &rhoL, const double &vL, const double &wL, const double &pL,
                                 const double &rhoR, const double &vR, const double &wR, const double &pR,
                                 Flux &fl, double &tempP, double &tempU)
    {
        if (tempU >= 0.0)
        { // pick left side
            if (tempP > pL)
                fl.updateFromPrimatives(rhoL * (((tempP / pL) + G6) / (G6 * (tempP / pL) + 1)),
                                        tempU, vL, wL, tempP);
            else
                fl.updateFromPrimatives(rhoL * pow((tempP / pL), 1 / G), tempU, vL, wL, tempP);
        }
        else
        { // pick right side
            if (tempP > pR)
                fl.updateFromPrimatives(rhoR * (((tempP / pR) + G6) / (G6 * (tempP / pR) + 1)),
                                        tempU, vR, wR, tempP);
            else
                fl.updateFromPrimatives(rhoR * pow((tempP / pR), 1 / G), tempU, vR, wR, tempP);
        }
        return true;
    }

    bool RiemannSolver::testVacuum(const double &rhoL, const double &uL, const double &vL, const double &wL,
                                   const double &aL, const double &pL, const double &rhoR, const double &uR,
                                   const double &vR, const double &wR, const double &aR, const double &pR,
                                   Flux &fl)
    {
        if ((rhoL == 0.0 || pL == 0.0) || (rhoR == 0.0 || pR == 0.0))
        {
            if (((rhoL == 0.0 || pL == 0.0) && (rhoR != 0.0 || pR != 0.0)) || (((rhoR != 0.0 || pR != 0.0) && (rhoL != 0.0 || pL != 0.0)) && 0 >= (uR - aR * G4))) // vacuum left not right
            {
                if (0 >= (uR + aR))
                {
                    std::cout << "W_R" << std::endl;
                    fl.updateFromPrimatives(rhoR, uR, vR, wR, pR);
                    // calca();
                }
                else if (0 <= (uR - aR * G4))
                {
                    std::cout << "W_L" << std::endl;
                    fl.updateFromPrimatives(rhoL, uL, vL, wL, pL);
                }
                else
                {
                    std::cout << "W_RFan" << std::endl;
                    fl.updateFromPrimatives(rhoR * pow(G5 - G6 * uR / aR, G4), G5 * (-aR + G7 * uR), vR * pow(G5 - G6 * uR / aR, G4), wR * pow(G5 - G6 * uR / aR, G4), pR * pow(G5 - G6 * uR / aR, G3));
                    // calca();
                }
                return true;
            }
            else if (((rhoR == 0.0 || pR == 0.0) && (rhoL != 0.0 || pL != 0.0)) || (((rhoR != 0.0 || pR != 0.0) && (rhoL != 0.0 || pL != 0.0)) && 0 <= (uL + aL * G4))) // vacuum right not left
            {
                if (0 <= (uL - aL))
                {
                    std::cout << "W_L" << std::endl;
                    fl.updateFromPrimatives(rhoL, uL, vL, wL, pL);
                    // calca();
                }
                else if (0 >= (uL + aL * G4))
                {
                    std::cout << "W_R" << std::endl;
                    fl.updateFromPrimatives(rhoR, uR, vR, wR, pR);
                }
                else
                {
                    std::cout << "W_LFan" << std::endl;
                    fl.updateFromPrimatives(rhoL * pow(G5 + G6 * uL / aL, G4), G5 * (aL + G7 * uL), vL * pow(G5 + G6 * uL / aL, G4), wL * pow(G5 + G6 * uL / aL, G4), pL * pow(G5 + G6 * uL / aL, G3));
                    // calca();
                }
                return true;
            }
            else if (((rhoR == 0.0 || pR == 0.0) && (rhoL == 0.0 || pL == 0.0)) || ((uL + aL * G4) < 0 && (uR - aR * G4) > 0)) // vacuum left and right
            {
                std::cout << "W_0" << std::endl;
                fl.updateFromPrimatives(0.0, 0.5 * (uL + uR), 0.5 * (vL + vR), 0.5 * (wL + wR), 0.0);
                return true;
            }
        }
        return false;
    }

    bool RiemannSolver::pickStartVal(const int errorStage, const double &rhoL, const double &uL,
                                     const double &aL, const double &pL, const double &rhoR, const double &uR,
                                     const double &aR, const double &pR, double &tempP)
    {
        // Implementation of pickStartVal...
        double p_PV = 0.5 * (pL + pR) + 0.5 * (uL - uR) * 0.25 * (rhoL + rhoR) * (aL + aR);
        p_PV = std::max(0.0, p_PV);

        switch (errorStage)
        {
        case 0:
        {
            double p_min = std::min(pL, pR);
            double p_max = std::max(pL, pR);
            double q_max = p_max / p_min;
            if (q_max < 2.0 && (p_min < p_PV && p_PV < p_max))
            {
                // Select PVRS Riemann solver
                tempP = p_PV;
            }
            else if (p_PV < p_min)
            {
                // Select Two-Rarefaction Riemann solver
                double p_q = pow(pL / pR, G1);
                double u_m = (p_q * uL / aL + uR / aR + G4 * (p_q - 1.0)) / (p_q / aL + 1 / aR);
                double p_TL = 1 + G7 * (uL - u_m) / aL;
                double p_TR = 1 + G7 * (u_m - uR) / aR;
                tempP = 0.5 * (pL * pow(p_TL, G3) + pR * pow(p_TR, G3));
            }
            else
            {
                // Select Two-Shock Riemann solver with
                // PVRS as estimate
                double ge_L = sqrt((G5 / rhoL) / (G6 * pL + p_PV));
                double ge_R = sqrt((G5 / rhoR) / (G6 * pR + p_PV));
                tempP = (ge_L * pL + ge_R * pR - (uR - uL)) / (ge_L + ge_R);
            }
            break;
        }
        case 1:
        {
            std::cout << "went to error stage 1" << std::endl;
            tempP = pow((aL + aR - 0.5 * G8 * (uR - uL)) / (aL / pow(pL, G1) + aR / pow(pR, G1)), G3);
            break;
        }
        case 2:
        {
            tempP = 0.5 * (pL + pR);
            break;
        }
        case 3:
        {
            double p_PV = 0.5 * (pL + pR) + 0.5 * (uL - uR) * 0.5 * (rhoL + rhoR) * 0.5 * (aL + aR);
            tempP = p_PV;
            break;
        }
        case 4:
        {
            double p_PV = 0.5 * (pL + pR) + 0.5 * (uL - uR) * 0.5 * (rhoL + rhoR) * 0.5 * (aL + aR);
            if (p_PV > TOL)
                tempP = p_PV;
            else
                tempP = TOL;
            double AL = 2 / ((G + 1) * rhoL);
            double BL = pL * (G - 1) / (G + 1);
            double AR = 2 / ((G + 1) * rhoR);
            double BR = pR * (G - 1) / (G + 1);
            double gL = pow(AL / (tempP + BL), 0.5);
            double gR = pow(AR / (tempP + BR), 0.5);
            double p_TS = (gL * pL + gR * pR - (uR - uL)) / (gL + gR);
            tempP = p_TS;
            break;
        }
        case 5:
        {
            tempP = p_PV;
            break;
        }
        case 6:
        {
            // Select Two-Rarefaction Riemann solver
            double p_q = pow(pL / pR, G1);
            double u_m = (p_q * uL / aL + uR / aR + G4 * (p_q - 1.0)) / (p_q / aL + 1 / aR);
            double p_TL = 1 + G7 * (uL - u_m) / aL;
            double p_TR = 1 + G7 * (u_m - uR) / aR;
            tempP = 0.5 * (pL * pow(p_TL, G3) + pR * pow(p_TR, G3));
            break;
        }
        case 7:
        {
            // Select Two-Shock Riemann solver with
            // PVRS as estimate
            double ge_L = sqrt((G5 / rhoL) / (G6 * pL + p_PV));
            double ge_R = sqrt((G5 / rhoR) / (G6 * pR + p_PV));
            tempP = (ge_L * pL + ge_R * pR - (uR - uL)) / (ge_L + ge_R);
            break;
        }
        case 8:
        {
            tempP = 1 / (rhoL * aL + rhoR * aR) * (rhoR * aR * pL + rhoL * aL * pR + rhoL * aL * rhoR * aR * (uL - aR));
            break;
        }
        case 9:
        {
            tempP = 1 * pow(10, -6);
            break;
        }
        default:
        {
            std::cout << "Error converging on p in x" << std::endl;
            std::cout << "pL: " << pL << ", rhoL: " << rhoL << ", uL: " << uL << ", aL: " << aL << ", pR: " << pR << ", rhoR: " << rhoR << ", uR: " << uR << ", aR: " << aR << std::endl;
            return false; // the timestep is probably too small
        }
        }
        return true;
    }

    bool RiemannSolver::iterateP(const double &rhoL, const double &uL, const double &aL, const double &pL,
                                 const double &rhoR, const double &uR, const double &aR, const double &pR,
                                 double &tempP, double &tempU)
    {
        bool iterate = true;
        double fsL, fsR, d_fsL, d_fsR, change;
        int errorStage = 0;

        while (iterate) // loop to try all the initial p values
        {
            iterate = pickStartVal(errorStage, rhoL, uL, aL, pL, rhoR, uR, aR, pR, tempP);
            if (!iterate)
                return false;
            int count = 0;
            // std::cout << "p guessed: " << p<< std::endl;
            while (iterate) // loop to iterate the p value
            {
                double A1 = 2 / ((G + 1) * rhoL);
                double B1 = pL * (G - 1) / (G + 1);
                double A2 = 2 / ((G + 1) * rhoR);
                double B2 = pR * (G - 1) / (G + 1);
                double p1 = A1 / (tempP + B1);
                double p2 = tempP / pL;
                double p3 = p2;
                double p4 = tempP / pR;
                double p5 = p4;
                double p6 = A2 / (tempP + B2);

                if (p1 < 0 || p2 < 0 || p3 < 0 || p4 < 0 || p5 < 0 || p6 < 0)
                {
                    errorStage++;
                    break; // exit the iterate loop to try next p
                }

                p1 = pow(p1, 0.5);
                p2 = pow(p2, G1);
                p3 = pow(p3, -G2);
                p4 = pow(p4, G1);
                p5 = pow(p5, -G2);
                p6 = pow(p6, 0.5);

                if (tempP > pL) // shock
                {
                    fsL = (tempP - pL) * p1;
                    d_fsL = p1 * (1 - (tempP - pL) / (2 * (B1 + tempP)));
                }
                else // expansion
                {
                    fsL = aL * G4 * (p2 - 1);
                    d_fsL = 1 / (rhoL * aL) * p3;
                }

                if (tempP > pR) // shock
                {
                    fsR = (tempP - pR) * p6;
                    d_fsR = p6 * (1 - (tempP - pR) / (2 * (B2 + tempP)));
                }
                else // expansion
                {
                    fsR = aR * G4 * (p4 - 1);
                    d_fsR = 1 / (rhoR * aR) * p5;
                }

                double f_ = fsL + fsR - uL + uR;
                double d_f = d_fsL + d_fsR;
                change = f_ / d_f;
                //    std::cout << f << ", " << d_f << ", " << change << std::endl;
                tempP = tempP - change; // Update new estimate of p*
                count++;

                if (TOL >= 2 * fabs(change / (change + 2 * tempP))) // iteration limit (slightly different to notes as abs of entire rhs)
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
        tempU = 0.5 * (uL + uR) + 0.5 * (fsR - fsL); // u*
        return true;
    }

} // namespace fluid