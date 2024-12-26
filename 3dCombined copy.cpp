#include <cmath>
#include <iostream>

const double gma = 1.4;
const double TOL = 0.000001;

struct values
{
    double rho, u, v, w, p, a;
};

struct boundary
{
    values l, r, f;
};

void calca(boundary &bD)
{
    bD.f.a = sqrt((gma * bD.f.p) / bD.f.rho);
}

bool pickSide(boundary &bD)
{
    // pick rho value depending on the side of the discontinuity
    if (bD.f.u >= 0.0) // pick left side
    {
        if (bD.f.p > bD.l.p)
            bD.f.rho = bD.l.rho * (((bD.f.p / bD.l.p) + ((gma - 1) / (gma + 1))) / (((gma - 1) / (gma + 1)) * (bD.f.p / bD.l.p) + 1));
        else
            bD.f.rho = bD.l.rho * pow((bD.f.p / bD.l.p), 1 / gma);
        bD.f.v = bD.l.v;
        bD.f.w = bD.l.w;
    }
    else // pick right side
    {
        if (bD.f.p > bD.r.p)
            bD.f.rho = bD.r.rho * (((bD.f.p / bD.r.p) + ((gma - 1) / (gma + 1))) / (((gma - 1) / (gma + 1)) * (bD.f.p / bD.r.p) + 1));
        else
            bD.f.rho = bD.r.rho * pow((bD.f.p / bD.r.p), 1 / gma);
        bD.f.v = bD.r.v;
        bD.f.w = bD.r.w;
    }
    if (std::isnan(bD.f.rho))
    {
        std::cout << "rhoerror u: " << bD.f.u << ", rhoL: " << bD.l.rho << ", rhoR: " << bD.r.rho << ", vL: " << bD.l.v << ", vR: " << bD.r.v << ", pL: " << bD.l.p << ", pR: " << bD.r.p << ", p:" << bD.f.p << std::endl;
        return false;
    }
    calca(bD);
    return true;
}

bool testVacuum(boundary &bD)
{
    if ((bD.l.rho == 0.0 || bD.l.p == 0.0) || (bD.r.rho == 0.0 || bD.r.p == 0.0))
    {
        if (((bD.l.rho == 0.0 || bD.l.p == 0.0) && (bD.r.rho != 0.0 || bD.r.p != 0.0)) || (((bD.r.rho != 0.0 || bD.r.p != 0.0) && (bD.l.rho != 0.0 || bD.l.p != 0.0)) && 0 >= (bD.r.u - 2 * bD.r.a / (gma - 1)))) // vacuum left not right
        {
            if (0 >= (bD.r.u + bD.r.a))
            {
                std::cout << "W_R" << std::endl;
                bD.f = bD.r;
                calca(bD);
            }
            else if (0 <= (bD.r.u - 2 * bD.r.a / (gma - 1)))
            {
                std::cout << "W_L" << std::endl;
                bD.f = bD.l;
            }
            else
            {
                std::cout << "W_RFan" << std::endl;
                bD.f.rho = bD.r.rho * pow(2 / (gma + 1) - (gma - 1) / (gma + 1) * (bD.r.u) / bD.r.a, 2 / (gma - 1));
                bD.f.u = 2 / (gma + 1) * (-bD.r.a + ((gma - 1) / 2) * bD.r.u);
                bD.f.v = bD.r.v * pow(2 / (gma + 1) - (gma - 1) / (gma + 1) * (bD.r.u) / bD.r.a, 2 / (gma - 1));
                bD.f.w = bD.r.w * pow(2 / (gma + 1) - (gma - 1) / (gma + 1) * (bD.r.u) / bD.r.a, 2 / (gma - 1));
                bD.f.p = bD.r.p * pow(2 / (gma + 1) - (gma - 1) / (gma + 1) * (bD.r.u) / bD.r.a, (2 * gma) / (gma - 1));
                calca(bD);
            }
            return true;
        }
        else if (((bD.r.rho == 0.0 || bD.r.p == 0.0) && (bD.l.rho != 0.0 || bD.l.p != 0.0)) || (((bD.r.rho != 0.0 || bD.r.p != 0.0) && (bD.l.rho != 0.0 || bD.l.p != 0.0)) && 0 <= (bD.l.u + 2 * bD.l.a / (gma - 1)))) // vacuum right not left
        {
            if (0 <= (bD.l.u - bD.l.a))
            {
                std::cout << "W_L" << std::endl;
                bD.f = bD.l;
                calca(bD);
            }
            else if (0 >= (bD.l.u + 2 * bD.l.a / (gma - 1)))
            {
                std::cout << "W_R" << std::endl;
                bD.f = bD.r;
            }
            else
            {
                std::cout << "W_LFan" << std::endl;
                bD.f.rho = bD.l.rho * pow(2 / (gma + 1) + (gma - 1) / (gma + 1) * (bD.l.u) / bD.l.a, 2 / (gma - 1));
                bD.f.u = 2 / (gma + 1) * (bD.l.a + ((gma - 1) / 2) * bD.l.u);
                bD.f.v = bD.l.v * pow(2 / (gma + 1) + (gma - 1) / (gma + 1) * (bD.l.u) / bD.l.a, 2 / (gma - 1));
                bD.f.w = bD.l.w * pow(2 / (gma + 1) + (gma - 1) / (gma + 1) * (bD.l.u) / bD.l.a, 2 / (gma - 1));
                bD.f.p = bD.l.p * pow(2 / (gma + 1) + (gma - 1) / (gma + 1) * (bD.l.u) / bD.l.a, (2 * gma) / (gma - 1));
                calca(bD);
            }
            return true;
        }
        else if (((bD.r.rho == 0.0 || bD.r.p == 0.0) && (bD.l.rho == 0.0 || bD.l.p == 0.0)) || ((bD.l.u + 2 * bD.l.a / (gma - 1)) < 0 && (bD.r.u - 2 * bD.r.a / (gma - 1)) > 0)) // vacuum left and right
        {
            std::cout << "W_0" << std::endl;
            bD.f.rho = 0.0;
            bD.f.p = 0.0;
            bD.f.u = 0.5 * (bD.l.u + bD.r.u); // estimate
            bD.f.v = 0.5 * (bD.l.v + bD.r.v); // estimate
            bD.f.w = 0.5 * (bD.l.w + bD.r.w); // estimate
            bD.f.a = 0.0;
            return true;
        }
    }
    return false;
}

bool pickStartVal(const int errorStage, boundary &bD)
{

    double G1 = (gma - 1) / (2 * gma);
    double G3 = 2 * gma / (gma - 1);
    double G4 = 2 / (gma - 1);
    double G5 = 2 / (gma + 1);
    double G6 = (gma - 1) / (gma + 1);
    double G7 = (gma - 1) / 2;

    double c_up = 0.25 * (bD.l.rho + bD.r.rho) * (bD.l.a + bD.r.a);
    double p_PV = 0.5 * (bD.l.p + bD.r.p) + 0.5 * (bD.l.u - bD.r.u) * c_up;
    p_PV = std::max(0.0, p_PV);

    // different guesses for p
    if (errorStage == 0)
    {
        double p_min = std::min(bD.l.p, bD.r.p);
        double p_max = std::max(bD.l.p, bD.r.p);
        double q_max = p_max / p_min;
        if (q_max < 2.0 && (p_min < p_PV && p_PV < p_max))
        {
            // Select PVRS Riemann solver
            bD.f.p = p_PV;
        }
        else if (p_PV < p_min)
        {
            // Select Two-Rarefaction Riemann solver
            double p_q = pow(bD.l.p / bD.r.p, G1);
            double u_m = (p_q * bD.l.u / bD.l.a + bD.r.u / bD.r.a + G4 * (p_q - 1.0)) / (p_q / bD.l.a + 1 / bD.r.a);
            double p_TL = 1 + G7 * (bD.l.u - u_m) / bD.l.a;
            double p_TR = 1 + G7 * (u_m - bD.r.u) / bD.r.a;
            bD.f.p = 0.5 * (bD.l.p * pow(p_TL, G3) + bD.r.p * pow(p_TR, G3));
        }
        else
        {
            // Select Two-Shock Riemann solver with
            // PVRS as estimate
            double ge_L = sqrt((G5 / bD.l.rho) / (G6 * bD.l.p + p_PV));
            double ge_R = sqrt((G5 / bD.r.rho) / (G6 * bD.r.p + p_PV));
            bD.f.p = (ge_L * bD.l.p + ge_R * bD.r.p - (bD.r.u - bD.l.u)) / (ge_L + ge_R);
        }
    }
    else if (errorStage == 1)
    {
        std::cout << "went to error stage 1" << std::endl;
        bD.f.p = pow((bD.l.a + bD.r.a - 0.5 * (gma - 1) * (bD.r.u - bD.l.u)) / (bD.l.a / pow(bD.l.p, (gma - 1) / (2 * gma)) + bD.r.a / pow(bD.r.p, (gma - 1) / (2 * gma))), (2 * gma) / (gma - 1));
    }
    else if (errorStage == 2)
    {
        bD.f.p = 0.5 * (bD.l.p + bD.r.p);
    }
    else if (errorStage == 3)
    {
        double p_PV = 0.5 * (bD.l.p + bD.r.p) + 0.5 * (bD.l.u - bD.r.u) * 0.5 * (bD.l.rho + bD.r.rho) * 0.5 * (bD.l.a + bD.r.a);
        // if (p_PV > TOL)
        bD.f.p = p_PV;
        //  else
        //     p = TOL;
    }
    else if (errorStage == 4)
    {
        double p_PV = 0.5 * (bD.l.p + bD.r.p) + 0.5 * (bD.l.u - bD.r.u) * 0.5 * (bD.l.rho + bD.r.rho) * 0.5 * (bD.l.a + bD.r.a);
        if (p_PV > TOL)
            bD.f.p = p_PV;
        else
            bD.f.p = TOL;
        double AL = 2 / ((gma + 1) * bD.l.rho);
        double BL = bD.l.p * (gma - 1) / (gma + 1);
        double AR = 2 / ((gma + 1) * bD.r.rho);
        double BR = bD.r.p * (gma - 1) / (gma + 1);
        double gL = pow(AL / (bD.f.p + BL), 0.5);
        double gR = pow(AR / (bD.f.p + BR), 0.5);
        double p_TS = (gL * bD.l.p + gR * bD.r.p - (bD.r.u - bD.l.u)) / (gL + gR);
        //   if (p_TS > TOL)
        bD.f.p = p_TS;
        //  else
        //      p = TOL;
    }
    else if (errorStage == 5)
    {
        bD.f.p = p_PV;
    }
    else if (errorStage == 6)
    {
        // Select Two-Rarefaction Riemann solver
        double p_q = pow(bD.l.p / bD.r.p, G1);
        double u_m = (p_q * bD.l.u / bD.l.a + bD.r.u / bD.r.a + G4 * (p_q - 1.0)) / (p_q / bD.l.a + 1 / bD.r.a);
        double p_TL = 1 + G7 * (bD.l.u - u_m) / bD.l.a;
        double p_TR = 1 + G7 * (u_m - bD.r.u) / bD.r.a;
        bD.f.p = 0.5 * (bD.l.p * pow(p_TL, G3) + bD.r.p * pow(p_TR, G3));
    }
    else if (errorStage == 7)
    {
        // Select Two-Shock Riemann solver with
        // PVRS as estimate
        double ge_L = sqrt((G5 / bD.l.rho) / (G6 * bD.l.p + p_PV));
        double ge_R = sqrt((G5 / bD.r.rho) / (G6 * bD.r.p + p_PV));
        bD.f.p = (ge_L * bD.l.p + ge_R * bD.r.p - (bD.r.u - bD.l.u)) / (ge_L + ge_R);
    }
    else if (errorStage == 8)
    {
        bD.f.p = 1 / (bD.l.rho * bD.l.a + bD.r.rho * bD.r.a) * (bD.r.rho * bD.r.a * bD.l.p + bD.l.rho * bD.l.a * bD.r.p + bD.l.rho * bD.l.a * bD.r.rho * bD.r.a * (bD.l.u - bD.r.a));
    }
    else if (errorStage == 9)
    {
        bD.f.p = 1 * pow(10, -6);
    }
    else
    {
        std::cout << "Error converging on p in x" << std::endl;
        std::cout << "pL: " << bD.l.p << ", rhoL: " << bD.l.rho << ", uL: " << bD.l.u << ", aL: " << bD.l.a << ", pR: " << bD.r.p << ", rhoR: " << bD.r.rho << ", uR: " << bD.r.u << ", aR: " << bD.r.a << std::endl;
        return false; // this doesnt actually stop the program but if you see that in the console the timestep is probably too small
    }
    return true;
}

bool iterateP(boundary &bD)
{
    bool iterate = true;
    double fsL, fsR, d_fsL, d_fsR, change;
    int errorStage = 0;

    while (iterate) // loop to try all the initial p values
    {
        iterate = pickStartVal(errorStage, bD);
        if (!iterate)
            return false;
        int count = 0;
        // std::cout << "p guessed: " << p<< std::endl;
        while (iterate) // loop to iterate the p value
        {
            double A1 = 2 / ((gma + 1) * bD.l.rho);
            double B1 = bD.l.p * (gma - 1) / (gma + 1);
            double A2 = 2 / ((gma + 1) * bD.r.rho);
            double B2 = bD.r.p * (gma - 1) / (gma + 1);
            double p1 = A1 / (bD.f.p + B1);
            double p2 = bD.f.p / bD.l.p;
            double p3 = p2;
            double p4 = bD.f.p / bD.r.p;
            double p5 = p4;
            double p6 = A2 / (bD.f.p + B2);

            if (p1 < 0 || p2 < 0 || p3 < 0 || p4 < 0 || p5 < 0 || p6 < 0)
            {
                errorStage++;
                break; // exit the iterate loop to try next p
            }

            p1 = pow(p1, 0.5);
            p2 = pow(p2, (gma - 1) / (2 * gma));
            p3 = pow(p3, -(gma + 1) / (2 * gma));
            p4 = pow(p4, (gma - 1) / (2 * gma));
            p5 = pow(p5, -(gma + 1) / (2 * gma));
            p6 = pow(p6, 0.5);

            if (bD.f.p > bD.l.p) // shock
            {
                fsL = (bD.f.p - bD.l.p) * p1;
                d_fsL = p1 * (1 - (bD.f.p - bD.l.p) / (2 * (B1 + bD.f.p)));
            }
            else // expansion
            {
                fsL = 2 * bD.l.a / (gma - 1) * (p2 - 1);
                d_fsL = 1 / (bD.l.rho * bD.l.a) * p3;
            }

            if (bD.f.p > bD.r.p) // shock
            {
                fsR = (bD.f.p - bD.r.p) * p6;
                d_fsR = p6 * (1 - (bD.f.p - bD.r.p) / (2 * (B2 + bD.f.p)));
            }
            else // expansion
            {
                fsR = 2 * bD.r.a / (gma - 1) * (p4 - 1);
                d_fsR = 1 / (bD.r.rho * bD.r.a) * p5;
            }

            double f = fsL + fsR - bD.l.u + bD.r.u;
            double d_f = d_fsL + d_fsR;
            change = f / d_f;
            //    std::cout << f << ", " << d_f << ", " << change << std::endl;
            bD.f.p = bD.f.p - change; // Update new estimate of p*
            count++;

            if (TOL >= 2 * fabs(change / (change + 2 * bD.f.p))) // iteration limit (slightly different to notes as abs of entire rhs)
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
    bD.f.u = 0.5 * (bD.l.u + bD.r.u) + 0.5 * (fsR - fsL); // u*
    return true;
}

bool findStar(boundary &bD) // find the values at the faces between 2 cells adjacent in x
{
    if (testVacuum(bD))
    {
        std::cout << "vacuum" << std::endl;
        return true; // a vacuum is generated, values found without iteration required NEED TO CHECK SPEED OF SOUND!!!!!
    }

    if (!iterateP(bD))
    {
        std::cout << "iteration failed" << std::endl;
        return false; // iteration failed, abort
    }

    if (!pickSide(bD))
    {
        std::cout << "pick side failed, rho is nan" << std::endl;
        return false; // rho is nan
    }
    return true;
}

int main()
{
    boundary bD;

    for (int caseNum = 1; caseNum < 7; caseNum++)
    {
        if (caseNum == 1)
        {
            bD.l.rho = 1.0;
            bD.l.u = 0.0;
            bD.l.v = 0.0;
            bD.l.w = 0.0;
            bD.l.p = 1.0;
            bD.r.rho = 0.125;
            bD.r.u = 0.0;
            bD.r.v = 0.0;
            bD.r.w = 0.0;
            bD.r.p = 0.1;
        }
        else if (caseNum == 2)
        {
            bD.l.p = 0.4;
            bD.l.rho = 1.0;
            bD.l.u = -2.0;
            bD.l.v = 0.0;
            bD.l.w = 0.0;
            bD.r.p = 0.4;
            bD.r.rho = 1.0;
            bD.r.u = 2.0;
            bD.r.v = 0.0;
            bD.r.w = 0.0;
        }
        else if (caseNum == 3)
        {
            bD.l.p = 1000.0;
            bD.l.rho = 1.0;
            bD.l.u = 0.0;
            bD.l.v = 0.0;
            bD.l.w = 0.0;
            bD.r.p = 0.01;
            bD.r.rho = 1.0;
            bD.r.u = 0.0;
            bD.r.v = 0.0;
            bD.r.w = 0.0;
        }
        else if (caseNum == 4)
        {
            bD.l.p = 0.01;
            bD.l.rho = 1.0;
            bD.l.u = 0.0;
            bD.l.v = 0.0;
            bD.l.w = 0.0;
            bD.r.p = 100.0;
            bD.r.rho = 1.0;
            bD.r.u = 0.0;
            bD.r.v = 0.0;
            bD.r.w = 0.0;
        }
        else if (caseNum == 5)
        {
            bD.l.p = 460.894;
            bD.l.rho = 5.99924;
            bD.l.u = 19.5975;
            bD.l.v = 0.0;
            bD.l.w = 0.0;
            bD.r.p = 46.0950;
            bD.r.rho = 5.99242;
            bD.r.u = -6.19633;
            bD.r.v = 0.0;
            bD.r.w = 0.0;
        }
        else if (caseNum == 6)
        {
            bD.l.p = 1.0;
            bD.l.rho = 1.0;
            bD.l.u = 0.0;
            bD.l.v = 0.0;
            bD.l.w = 0.0;
            bD.r.p = 0.0;
            bD.r.rho = 0.0;
            bD.r.u = 0.0;
            bD.r.v = 0.0;
            bD.r.w = 0.0;
        }

        if (bD.l.rho == 0.0 || bD.l.p == 0.0)
        {
            bD.l.a = 0.0;
        }
        else
        {
            bD.l.a = sqrt((gma * bD.l.p) / bD.l.rho);
        }
        if (bD.r.rho == 0.0 || bD.r.p == 0.0)
        {
            bD.r.a = 0.0;
        }
        else
        {
            bD.r.a = sqrt((gma * bD.r.p) / bD.r.rho);
        }
        findStar(bD);
        std::cout << "caseNum: " << caseNum << ", p: " << bD.f.p << ", u: " << bD.f.u << ", rho: " << bD.f.rho << std::endl;
        std::cout << "caseNum: " << caseNum << ", f1: " << bD.f.rho * bD.f.u << ", f2: " << bD.f.rho * bD.f.u * bD.f.u + bD.f.p << ", f3: " << bD.f.rho * bD.f.v * bD.f.u << ", f4: " << bD.f.rho * bD.f.w * bD.f.u << ", f5: " << bD.f.u * (bD.f.rho * (0.5 * (bD.f.u * bD.f.u + bD.f.v * bD.f.v + bD.f.w * bD.f.w) + bD.f.p / ((gma - 1) * bD.f.rho) + bD.f.p)) << std::endl;
    }
    return 0;
}