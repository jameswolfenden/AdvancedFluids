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
    a = sqrt((gamma * p) / rho);
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
    return rho * (0.5 * (u * u + v * v) + p / ((gamma - 1) * rho));
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
    return u * (rho * (0.5 * (u * u + v * v) + p / ((gamma - 1) * rho)) + p);
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
    return v * (rho * (0.5 * (u * u + v * v) + p / ((gamma - 1) * rho)) + p);
}

void Cell::updateConservatives(double u1, double u2, double u3, double u4) // update the primatives from new conservative values
{
    rho = u1;
    u = u2 / u1;
    v = u3 / u1;
    p = (gamma - 1) * (u4 - 0.5 * ((u2 * u2 + u3 * u3) / u1)); // these all need checking i guessed it ngl
    aCalc();
}
void Cell::updatePrimatives(double p, double rho, double u, double v)
{
    this->p = p;
    this->rho = rho;
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

    while (iterate) // loop to try all the initial p values
    {
        // different guesses for p
        if (errorStage == 0)
        {
            p = pow((sides[0]->a + sides[1]->a - 0.5 * (gamma - 1) * (sides[1]->u - sides[0]->u)) / (sides[0]->a / pow(sides[0]->p, (gamma - 1) / (2 * gamma)) + sides[1]->a / pow(sides[1]->p, (gamma - 1) / (2 * gamma))), (2 * gamma) / (gamma - 1));
        }
        else if (errorStage == 1)
        {
            p = 0.5 * (sides[0]->p + sides[1]->p);
        }
        else if (errorStage == 2)
        {
            double p_PV = 0.5 * (sides[0]->p + sides[1]->p) + 0.5 * (sides[0]->u - sides[1]->u) * 0.5 * (sides[0]->rho + sides[1]->rho) * 0.5 * (sides[0]->a + sides[1]->a);
            if (p_PV > TOL)
                p = p_PV;
            else
                p = TOL;
        }
        else
        {
            std::cout << "Error converging on p in x" << std::endl;
            return; // this doesnt actually stop the program but if you see that in the console the timestep is probably too small
        }
        count = 0;
        while (iterate) // loop to iterate the p value
        {
            errno = 0;
            for (int side = 0; side < 2; side++)
            {
                double A = 2 / ((gamma + 1) * sides[side]->rho);
                double B = sides[side]->p * (gamma - 1) / (gamma + 1);
                if (p > sides[side]->p) // shock
                {
                    fs[side] = (p - sides[side]->p) * pow(A / (p + B), 0.5);
                    d_fs[side] = pow(A / (p + B), 0.5) * (1 - (p - sides[side]->p) / (2 * (B + p)));
                }
                else // expansion
                {
                    fs[side] = 2 * sides[side]->a / (gamma - 1) * (pow(p / sides[side]->p, (gamma - 1) / (2 * gamma)) - 1);
                    d_fs[side] = 1 / (sides[side]->p * sides[side]->a) * pow(p / sides[side]->p, -(gamma + 1) / (2 * gamma));
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
        rho = sides[0]->rho * (((p / sides[0]->p) + ((gamma - 1) / (gamma + 1))) / (((gamma - 1) / (gamma + 1)) * (p / sides[0]->p) + 1));
        v = sides[0]->v;
    }
    else
    {
        rho = sides[1]->rho * (((p / sides[1]->p) + ((gamma - 1) / (gamma + 1))) / (((gamma - 1) / (gamma + 1)) * (p / sides[1]->p) + 1));
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

    while (iterate)
    {
        if (errorStage == 0)
        {
            p = pow((sides[0]->a + sides[1]->a - 0.5 * (gamma - 1) * (sides[1]->v - sides[0]->v)) / (sides[0]->a / pow(sides[0]->p, (gamma - 1) / (2 * gamma)) + sides[1]->a / pow(sides[1]->p, (gamma - 1) / (2 * gamma))), (2 * gamma) / (gamma - 1));
        }
        else if (errorStage == 1)
        {
            p = 0.5 * (sides[0]->p + sides[1]->p);
        }
        else if (errorStage == 2)
        {
            double p_PV = 0.5 * (sides[0]->p + sides[1]->p) + 0.5 * (sides[0]->v - sides[1]->v) * 0.5 * (sides[0]->rho + sides[1]->rho) * 0.5 * (sides[0]->a + sides[1]->a);
            if (p_PV > TOL)
                p = p_PV;
            else
                p = TOL;
        }
        else
        {
            std::cout << "Error converging on p in y" << std::endl;
            return;
        }

        count = 0;
        while (iterate)
        {
            errno = 0;
            for (int side = 0; side < 2; side++)
            {
                double A = 2 / ((gamma + 1) * sides[side]->rho);
                double B = sides[side]->p * (gamma - 1) / (gamma + 1);
                if (p > sides[side]->p) // shock
                {
                    fs[side] = (p - sides[side]->p) * pow(A / (p + B), 0.5);
                    d_fs[side] = pow(A / (p + B), 0.5) * (1 - (p - sides[side]->p) / (2 * (B + p)));
                }
                else // expansion
                {
                    fs[side] = 2 * sides[side]->a / (gamma - 1) * (pow(p / sides[side]->p, (gamma - 1) / (2 * gamma)) - 1);
                    d_fs[side] = 1 / (sides[side]->p * sides[side]->a) * pow(p / sides[side]->p, -(gamma + 1) / (2 * gamma));
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
        rho = sides[0]->rho * (((p / sides[0]->p) + ((gamma - 1) / (gamma + 1))) / (((gamma - 1) / (gamma + 1)) * (p / sides[0]->p) + 1));
        u = sides[0]->u;
    }
    else
    {
        rho = sides[1]->rho * (((p / sides[1]->p) + ((gamma - 1) / (gamma + 1))) / (((gamma - 1) / (gamma + 1)) * (p / sides[1]->p) + 1));
        u = sides[1]->u;
    }
    aCalc();
    return;
}
