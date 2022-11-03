#include "Cell.h"
#include <cmath>
#include "../fluidConsts.h"
#include <iostream>

Cell::Cell(double p, double rho, double u, double v)
{
    updatePrimatives(p, rho, u, v);
}
double Cell::aCalc() // find speed of sound
{
    if (rho == 0)
        a = 0;
    else
        a = sqrt((gammma * p) / rho);
    return a;
}
// the following are to get the conservatives and fluxes from the primatives
double Cell::u1()
{
    return rho;
}
double Cell::u2()
{
    return rho * u;
}
double Cell::u3()
{
    return rho * v;
}
double Cell::u4()
{
    if (rho == 0)
        return 0;
    else
        return rho * (0.5 * (u * u + v * v) + p / ((gammma - 1) * rho));
}
double Cell::f1()
{
    return rho * u;
}
double Cell::f2()
{
    return rho * u * u + p;
}
double Cell::f3()
{
    return rho * v * u;
}
double Cell::f4()
{
    if (rho == 0)
        return 0;
    else
        return u * (rho * (0.5 * (u * u + v * v) + p / ((gammma - 1) * rho)) + p);
}
double Cell::g1()
{
    return rho * v;
}
double Cell::g2()
{
    return rho * v * u;
}
double Cell::g3()
{
    return rho * v * v + p;
}
double Cell::g4()
{
    if (rho == 0)
        return 0;
    else
        return v * (rho * (0.5 * (u * u + v * v) + p / ((gammma - 1) * rho)) + p);
}

void Cell::updateConservatives(double u1, double u2, double u3, double u4) // update the primatives from new conservative values
{
    if (u1 == 0)
    {
        updatePrimatives(0, 0, 0, 0);
    }
    else
    {
        double rho_ = u1;
        double u_ = u2 / u1;
        double v_ = u3 / u1;
        double p_ = (gammma - 1) * (u4 - 0.5 * ((u2 * u2 + u3 * u3) / u1)); // these all need checking i guessed it ngl
        updatePrimatives(p_, rho_, u_, v_);
    }
}
void Cell::updatePrimatives(double p, double rho, double u, double v)
{
    if (p <= 0.0 || rho <= 0.0)
    {
        this->p = 0.0;
        this->rho = 0.0;
    }
    else
    {
        this->p = p;
        this->rho = rho;
    }

    this->u = u;
    this->v = v;
    aCalc();
}
void Cell::xFindStar(Cell *sides[]) // find the values at the faces between 2 cells adjacent in x
{

    double change;
    int count = 0;
    int errorStage = 0;
    bool iterate = true;
    double fs[2];
    double d_fs[2];

    if ((sides[0]->rho == 0.0 && sides[1]->rho != 0.0) || ((sides[1]->rho == 0.0 && sides[0]->rho == 0.0) && 0 >= (sides[1]->u - 2 * sides[1]->a / (gammma - 1)))) // vacuum left not right
    {
        if (0 >= (sides[1]->u + sides[1]->a))
        {
            rho = sides[1]->rho;
            u = sides[1]->u;
            v = sides[1]->v; // made up
            p = sides[1]->p;
        }
        else if (0 <= (sides[1]->u - 2 * sides[1]->a / (gammma - 1)))
        {
            rho = 0;
            p = 0;
            u = 0.5 * (sides[0]->u + sides[1]->u); // SOURCE: I MADE IT UP NEED TO ASK
            v = 0.5 * (sides[0]->v + sides[1]->v); // SOURCE: I MADE IT UP NEED TO ASK
        }
        else
        {
            rho = sides[1]->rho * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[1]->u) / sides[1]->a, 2 / (gammma - 1));
            u = 2 / (gammma + 1) * (-sides[1]->a + ((gammma - 1) / 2) * sides[1]->u);
            v = 2 / (gammma + 1) * (-sides[1]->a + ((gammma - 1) / 2) * sides[1]->v); // made up
            p = sides[1]->p * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[1]->u) / sides[1]->a, (2 * gammma) / (gammma - 1));
        }
        return;
    }
    else if ((sides[1]->rho == 0.0 && sides[0]->rho != 0.0) || ((sides[1]->rho == 0.0 && sides[0]->rho == 0.0) && 0 <= (sides[0]->u + 2 * sides[0]->a / (gammma - 1)))) // vacuum right not left
    {
        if (0 <= (sides[0]->u - sides[0]->a))
        {
            rho = sides[1]->rho;
            u = sides[1]->u;
            v = sides[1]->v; // made up
            p = sides[1]->p;
        }
        else if (0 >= (sides[0]->u + 2 * sides[0]->a / (gammma - 1)))
        {
            rho = 0;
            p = 0;
            u = 0.5 * (sides[0]->u + sides[1]->u); // SOURCE: I MADE IT UP NEED TO ASK
            v = 0.5 * (sides[0]->v + sides[1]->v); // SOURCE: I MADE IT UP NEED TO ASK
        }
        else
        {
            rho = sides[0]->rho * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[0]->u) / sides[0]->a, 2 / (gammma - 1));
            u = 2 / (gammma + 1) * (-sides[0]->a + ((gammma - 1) / 2) * sides[0]->u);
            v = 2 / (gammma + 1) * (-sides[0]->a + ((gammma - 1) / 2) * sides[0]->v); // made up
            p = sides[0]->p * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[0]->u) / sides[0]->a, (2 * gammma) / (gammma - 1));
        }
        return;
    }
    else if (sides[1]->rho == 0.0 && sides[0]->rho == 0.0) // vacuum left and right
    {
        rho = 0;
        p = 0;
        u = 0.5 * (sides[0]->u + sides[1]->u); // SOURCE: I MADE IT UP NEED TO ASK
        v = 0.5 * (sides[0]->v + sides[1]->v); // SOURCE: I MADE IT UP NEED TO ASK
        return;
    }

    while (iterate) // loop to try all the initial p values
    {
        // different guesses for p
        if (errorStage == 0)
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

            double c_up = 0.25 * (sides[0]->rho + sides[1]->rho) * (sides[0]->a + sides[1]->a);
            double p_PV = 0.5 * (sides[0]->p + sides[1]->p) + 0.5 * (sides[0]->u - sides[1]->u) * c_up;
            p_PV = std::max(0.0, p_PV);
            double p_min = std::min(sides[0]->p, sides[1]->p);
            double p_max = std::max(sides[0]->p, sides[1]->p);
            double q_max = p_max / p_min;

            if (q_max < q_user && (p_min < p_PV && p_PV < p_max))
            {
                // Select PVRS Riemann solver
                p = p_PV;
            }
            else if (p_PV < p_min)
            {
                // Select Two-Rarefaction Riemann solver
                double p_q = pow(sides[0]->p / sides[1]->p, G1);
                double u_m = (p_q * sides[0]->u / sides[0]->a + sides[1]->u / sides[1]->a + G4 * (p_q - 1.0)) / (p_q / sides[0]->a + 1 / sides[1]->a);
                double p_TL = 1 + G7 * (sides[0]->u - u_m) / sides[0]->a;
                double p_TR = 1 + G7 * (u_m - sides[1]->u) / sides[1]->a;
                p = 0.5 * (sides[0]->p * pow(p_TL, G3) + sides[1]->p * pow(p_TR, G3));
            }
            else
            {
                // Select Two-Shock Riemann solver with
                // PVRS as estimate
                double ge_L = sqrt((G5 / sides[0]->rho) / (G6 * sides[0]->p + p_PV));
                double ge_R = sqrt((G5 / sides[1]->rho) / (G6 * sides[1]->p + p_PV));
                p = (ge_L * sides[0]->p + ge_R * sides[1]->p - (sides[1]->u - sides[0]->u)) / (ge_L + ge_R);
            }
        }
        else if (errorStage == 1)
        {
            p = pow((sides[0]->a + sides[1]->a - 0.5 * (gammma - 1) * (sides[1]->u - sides[0]->u)) / (sides[0]->a / pow(sides[0]->p, (gammma - 1) / (2 * gammma)) + sides[1]->a / pow(sides[1]->p, (gammma - 1) / (2 * gammma))), (2 * gammma) / (gammma - 1));
        }
        else if (errorStage == 2)
        {
            p = 0.5 * (sides[0]->p + sides[1]->p);
        }
        else if (errorStage == 3)
        {
            double p_PV = 0.5 * (sides[0]->p + sides[1]->p) + 0.5 * (sides[0]->u - sides[1]->u) * 0.5 * (sides[0]->rho + sides[1]->rho) * 0.5 * (sides[0]->a + sides[1]->a);
            if (p_PV > TOL)
                p = p_PV;
            else
                p = TOL;
        }
        else if (errorStage == 4)
        {
            double p_PV = 0.5 * (sides[0]->p + sides[1]->p) + 0.5 * (sides[0]->u - sides[1]->u) * 0.5 * (sides[0]->rho + sides[1]->rho) * 0.5 * (sides[0]->a + sides[1]->a);
            if (p_PV > TOL)
                p = p_PV;
            else
                p = TOL;
            double AL = 2 / ((gammma + 1) * sides[0]->rho);
            double BL = sides[0]->p * (gammma - 1) / (gammma + 1);
            double AR = 2 / ((gammma + 1) * sides[1]->rho);
            double BR = sides[1]->p * (gammma - 1) / (gammma + 1);
            double gL = pow(AL / (p + BL), 0.5);
            double gR = pow(AR / (p + BR), 0.5);
            double p_TS = (gL * sides[0]->p + gR * sides[1]->p - (sides[1]->u - sides[0]->u)) / (gL + gR);
            if (p_TS > TOL)
                p = p_TS;
            else
                p = TOL;
        }
        else
        {
            std::cout << "Error converging on p in x" << std::endl;
            std::cout << "pL: " << sides[0]->p << ", rhoL: " << sides[0]->rho << ", uL: " << sides[0]->u << ", pR: " << sides[1]->p << ", rhoR: " << sides[1]->rho << ", uR: " << sides[1]->u << std::endl;
            return; // this doesnt actually stop the program but if you see that in the console the timestep is probably too small
        }
        count = 0;
        while (iterate) // loop to iterate the p value
        {
            errno = 0;
            for (int side = 0; side < 2; side++)
            {
                double A = 2 / ((gammma + 1) * sides[side]->rho);
                double B = sides[side]->p * (gammma - 1) / (gammma + 1);
                if (p > sides[side]->p) // shock
                {
                    fs[side] = (p - sides[side]->p) * pow(A / (p + B), 0.5);
                    d_fs[side] = pow(A / (p + B), 0.5) * (1 - (p - sides[side]->p) / (2 * (B + p)));
                }
                else // expansion
                {
                    fs[side] = 2 * sides[side]->a / (gammma - 1) * (pow(p / sides[side]->p, (gammma - 1) / (2 * gammma)) - 1);
                    d_fs[side] = 1 / (sides[side]->rho * sides[side]->a) * pow(p / sides[side]->p, -(gammma + 1) / (2 * gammma));
                }
            }
            if (errno != 0) // p is probably negative giving an issue with pow, we need the next guess of p
            {
                errorStage++;
                break; // exit the iterate loop to try next p
            }

            double f = fs[0] + fs[1] - sides[0]->u + sides[1]->u;
            double d_f = d_fs[0] + d_fs[1];
            change = f / d_f;
            p = p - change; // Update new estimate of p*
            count++;

            if (TOL >= 2 * fabs(change / (change + 2 * p))) // iteration limit (slightly different to notes as abs of entire rhs)
            {
                iterate = false; // we have converged
            }
        }
    }
    u = 0.5 * (sides[0]->u + sides[1]->u) + 0.5 * (fs[1] - fs[0]); // u*
    if (u >= 0)                                                    // pick rho value depending on the side of the discontinuity
    {
        if (p > sides[0]->p)
        {
            rho = sides[0]->rho * (((p / sides[0]->p) + ((gammma - 1) / (gammma + 1))) / (((gammma - 1) / (gammma + 1)) * (p / sides[0]->p) + 1));
        }
        else
        {
            rho = sides[0]->rho * pow((p / sides[0]->p), 1 / gammma);
        }
        v = sides[0]->v;
    }
    else
    {
        if (p > sides[1]->p)
        {
            rho = sides[1]->rho * (((p / sides[1]->p) + ((gammma - 1) / (gammma + 1))) / (((gammma - 1) / (gammma + 1)) * (p / sides[1]->p) + 1));
        }
        else
        {
            rho = sides[1]->rho * pow((p / sides[1]->p), 1 / gammma);
        }
        v = sides[1]->v;
    }
    aCalc();
    return;
}

void Cell::yFindStar(Cell *sides[]) // same as the x but with x and y and u and v swapped, could probably combine somehow but it works, see comments on x
{

    double change;
    int count = 0;
    int errorStage = 0;
    bool iterate = true;
    double fs[2];
    double d_fs[2];

    if ((sides[0]->rho == 0.0 && sides[1]->rho != 0.0) || ((sides[1]->rho == 0.0 && sides[0]->rho == 0.0) && 0 >= (sides[1]->v - 2 * sides[1]->a / (gammma - 1)))) // vacuum left not right
    {
        if (0 >= (sides[1]->v + sides[1]->a))
        {
            rho = sides[1]->rho;
            v = sides[1]->v;
            u = sides[1]->u; // made up
            p = sides[1]->p;
        }
        else if (0 <= (sides[1]->u - 2 * sides[1]->a / (gammma - 1)))
        {
            rho = 0;
            p = 0;
            v = 0.5 * (sides[0]->v + sides[1]->v); // SOURCE: I MADE IT UP NEED TO ASK
            u = 0.5 * (sides[0]->u + sides[1]->u); // SOURCE: I MADE IT UP NEED TO ASK
        }
        else
        {
            rho = sides[1]->rho * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[1]->v) / sides[1]->a, 2 / (gammma - 1));
            v = 2 / (gammma + 1) * (-sides[1]->a + ((gammma - 1) / 2) * sides[1]->v); // made up
            u = 2 / (gammma + 1) * (-sides[1]->a + ((gammma - 1) / 2) * sides[1]->u);
            p = sides[1]->p * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[1]->v) / sides[1]->a, (2 * gammma) / (gammma - 1));
        }
        return;
    }
    else if ((sides[1]->rho == 0.0 && sides[0]->rho != 0.0) || ((sides[1]->rho == 0.0 && sides[0]->rho == 0.0) && 0 <= (sides[0]->v + 2 * sides[0]->a / (gammma - 1)))) // vacuum right not left
    {
        if (0 <= (sides[0]->v - sides[0]->a))
        {
            rho = sides[1]->rho;
            v = sides[1]->v; // made up
            u = sides[1]->u;
            p = sides[1]->p;
        }
        else if (0 >= (sides[0]->v + 2 * sides[0]->a / (gammma - 1)))
        {
            rho = 0;
            p = 0;
            v = 0.5 * (sides[0]->v + sides[1]->v); // SOURCE: I MADE IT UP NEED TO ASK
            u = 0.5 * (sides[0]->u + sides[1]->u); // SOURCE: I MADE IT UP NEED TO ASK
        }
        else
        {
            rho = sides[0]->rho * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[0]->v) / sides[0]->a, 2 / (gammma - 1));
            v = 2 / (gammma + 1) * (-sides[0]->a + ((gammma - 1) / 2) * sides[0]->v); // made up
            u = 2 / (gammma + 1) * (-sides[0]->a + ((gammma - 1) / 2) * sides[0]->u);
            p = sides[0]->p * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[0]->v) / sides[0]->a, (2 * gammma) / (gammma - 1));
        }
        return;
    }
    else if (sides[1]->rho == 0.0 && sides[0]->rho == 0.0) // vacuum left and right
    {
        rho = 0;
        p = 0;
        v = 0.5 * (sides[0]->v + sides[1]->v); // SOURCE: I MADE IT UP NEED TO ASK
        u = 0.5 * (sides[0]->u + sides[1]->u); // SOURCE: I MADE IT UP NEED TO ASK
        return;
    }

    while (iterate)
    {
        if (errorStage == 0)
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

            double c_up = 0.25 * (sides[0]->rho + sides[1]->rho) * (sides[0]->a + sides[1]->a);
            double p_PV = 0.5 * (sides[0]->p + sides[1]->p) + 0.5 * (sides[0]->v - sides[1]->v) * c_up;
            p_PV = std::max(0.0, p_PV);
            double p_min = std::min(sides[0]->p, sides[1]->p);
            double p_max = std::max(sides[0]->p, sides[1]->p);
            double q_max = p_max / p_min;

            if (q_max < q_user && (p_min < p_PV && p_PV < p_max))
            {
                // Select PVRS Riemann solver
                p = p_PV;
            }
            else if (p_PV < p_min)
            {
                // Select Two-Rarefaction Riemann solver
                double p_q = pow(sides[0]->p / sides[1]->p, G1);
                double u_m = (p_q * sides[0]->v / sides[0]->a + sides[1]->v / sides[1]->a + G4 * (p_q - 1.0)) / (p_q / sides[0]->a + 1 / sides[1]->a);
                double p_TL = 1 + G7 * (sides[0]->v - u_m) / sides[0]->a;
                double p_TR = 1 + G7 * (u_m - sides[1]->v) / sides[1]->a;
                p = 0.5 * (sides[0]->p * pow(p_TL, G3) + sides[1]->p * pow(p_TR, G3));
            }
            else
            {
                // Select Two-Shock Riemann solver with
                // PVRS as estimate
                double ge_L = sqrt((G5 / sides[0]->rho) / (G6 * sides[0]->p + p_PV));
                double ge_R = sqrt((G5 / sides[1]->rho) / (G6 * sides[1]->p + p_PV));
                p = (ge_L * sides[0]->p + ge_R * sides[1]->p - (sides[1]->v - sides[0]->v)) / (ge_L + ge_R);
            }
        }
        else if (errorStage == 1)
        {
            p = pow((sides[0]->a + sides[1]->a - 0.5 * (gammma - 1) * (sides[1]->v - sides[0]->v)) / (sides[0]->a / pow(sides[0]->p, (gammma - 1) / (2 * gammma)) + sides[1]->a / pow(sides[1]->p, (gammma - 1) / (2 * gammma))), (2 * gammma) / (gammma - 1));
        }
        else if (errorStage == 2)
        {
            p = 0.5 * (sides[0]->p + sides[1]->p);
        }
        else if (errorStage == 3)
        {
            double p_PV = 0.5 * (sides[0]->p + sides[1]->p) + 0.5 * (sides[0]->v - sides[1]->v) * 0.5 * (sides[0]->rho + sides[1]->rho) * 0.5 * (sides[0]->a + sides[1]->a);
            if (p_PV > TOL)
                p = p_PV;
            else
                p = TOL;
        }
        else if (errorStage == 4)
        {
            double p_PV = 0.5 * (sides[0]->p + sides[1]->p) + 0.5 * (sides[0]->v - sides[1]->v) * 0.5 * (sides[0]->rho + sides[1]->rho) * 0.5 * (sides[0]->a + sides[1]->a);
            if (p_PV > TOL)
                p = p_PV;
            else
                p = TOL;
            double AL = 2 / ((gammma + 1) * sides[0]->rho);
            double BL = sides[0]->p * (gammma - 1) / (gammma + 1);
            double AR = 2 / ((gammma + 1) * sides[1]->rho);
            double BR = sides[1]->p * (gammma - 1) / (gammma + 1);
            double gL = pow(AL / (p + BL), 0.5);
            double gR = pow(AR / (p + BR), 0.5);
            double p_TS = (gL * sides[0]->p + gR * sides[1]->p - (sides[1]->v - sides[0]->v)) / (gL + gR);
            if (p_TS > TOL)
                p = p_TS;
            else
                p = TOL;
        }
        else
        {
            std::cout << "Error converging on p in y" << std::endl;
            std::cout << "pL: " << sides[0]->p << ", rhoL: " << sides[0]->rho << ", vL: " << sides[0]->v << ", pR: " << sides[1]->p << ", rhoR: " << sides[1]->rho << ", vR: " << sides[1]->v << std::endl;
            return;
        }

        count = 0;
        while (iterate)
        {
            errno = 0;
            for (int side = 0; side < 2; side++)
            {
                double A = 2 / ((gammma + 1) * sides[side]->rho);
                double B = sides[side]->p * (gammma - 1) / (gammma + 1);
                if (p > sides[side]->p) // shock
                {
                    fs[side] = (p - sides[side]->p) * pow(A / (p + B), 0.5);
                    d_fs[side] = pow(A / (p + B), 0.5) * (1 - (p - sides[side]->p) / (2 * (B + p)));
                }
                else // expansion
                {
                    fs[side] = 2 * sides[side]->a / (gammma - 1) * (pow(p / sides[side]->p, (gammma - 1) / (2 * gammma)) - 1);
                    d_fs[side] = 1 / (sides[side]->rho * sides[side]->a) * pow(p / sides[side]->p, -(gammma + 1) / (2 * gammma));
                }
            }
            if (errno != 0)
            {
                errorStage++;
                break;
            }
            double f = fs[0] + fs[1] - sides[0]->v + sides[1]->v;
            double d_f = d_fs[0] + d_fs[1];
            change = f / d_f;
            p = p - change; // Update new estimate of p*
            count++;
            if (TOL >= 2 * fabs(change / (change + 2 * p))) // iteration limit (slightly different to notes as abs of entire rhs)
            {
                iterate = false;
            }
        }
    }
    v = 0.5 * (sides[0]->v + sides[1]->v) + 0.5 * (fs[1] - fs[0]); // u*
    if (v >= 0)                                                    // pick rho value depending on the side of the discontinuity
    {
        if (p > sides[0]->p)
        {
            rho = sides[0]->rho * (((p / sides[0]->p) + ((gammma - 1) / (gammma + 1))) / (((gammma - 1) / (gammma + 1)) * (p / sides[0]->p) + 1));
        }
        else
        {
            rho = sides[0]->rho * pow((p / sides[0]->p), 1 / gammma);
        }
        u = sides[0]->u;
    }
    else
    {
        if (p > sides[1]->p)
        {
            rho = sides[1]->rho * (((p / sides[1]->p) + ((gammma - 1) / (gammma + 1))) / (((gammma - 1) / (gammma + 1)) * (p / sides[1]->p) + 1));
        }
        else
        {
            rho = sides[1]->rho * pow((p / sides[1]->p), 1 / gammma);
        }
        u = sides[1]->u;
    }
    aCalc();
    return;
}
