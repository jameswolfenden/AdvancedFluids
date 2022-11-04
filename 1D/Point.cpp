#include "Point.h"
#include <cmath>
#include "../fluidConsts.h"
#include <iostream>

Point::Point(double p, double rho, double u)
{
    updatePrimatives(p, rho, u);
}
double Point::aCalc()
{
    a = sqrt((gammma * p) / rho);
    return a;
}
double Point::u1()
{
    return rho;
}
double Point::u2()
{
    return rho * u;
}
double Point::u3()
{
    return rho * (0.5 * u * u + p / ((gammma - 1) * rho));
}
double Point::f1()
{
    return rho * u;
}
double Point::f2()
{
    return rho * u * u + p;
}
double Point::f3()
{
    return u * (rho * (0.5 * u * u + p / ((gammma - 1) * rho)) + p);
}
void Point::updateConservatives(double u1, double u2, double u3)
{
    rho = u1;
    u = u2 / u1;
    p = (gammma - 1) * (u3 - 0.5 * ((u2 * u2) / u1));
    aCalc();
}
void Point::updatePrimatives(double p, double rho, double u)
{
    this->p = p;
    this->rho = rho;
    this->u = u;
    aCalc();
}
void Point::findStar(Point *sides[])
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
            p = sides[1]->p;
        }
        else if (0 <= (sides[1]->u - 2 * sides[1]->a / (gammma - 1)))
        {
            rho = 0;
            p = 0;
            u = 0.5 * (sides[0]->u + sides[1]->u); // SOURCE: I MADE IT UP NEED TO ASK
        }
        else
        {
            rho = sides[1]->rho * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[1]->u) / sides[1]->a, 2 / (gammma - 1));
            u = 2 / (gammma + 1) * (-sides[1]->a + ((gammma - 1) / 2) * sides[1]->u);
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
            p = sides[1]->p;
        }
        else if (0 >= (sides[0]->u + 2 * sides[0]->a / (gammma - 1)))
        {
            rho = 0;
            p = 0;
            u = 0.5 * (sides[0]->u + sides[1]->u); // SOURCE: I MADE IT UP NEED TO ASK
        }
        else
        {
            rho = sides[0]->rho * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[0]->u) / sides[0]->a, 2 / (gammma - 1));
            u = 2 / (gammma + 1) * (-sides[0]->a + ((gammma - 1) / 2) * sides[0]->u);
            p = sides[0]->p * pow(2 / (gammma + 1) - (gammma - 1) / (gammma + 1) * (sides[0]->u) / sides[0]->a, (2 * gammma) / (gammma - 1));
        }
        return;
    }
    else if (sides[1]->rho == 0.0 && sides[0]->rho == 0.0) // vacuum left and right
    {
        rho = 0;
        p = 0;
        u = 0.5 * (sides[0]->u + sides[1]->u); // SOURCE: I MADE IT UP NEED TO ASK
        return;
    }

    while (iterate)
    {
        if (waveOutput)
            std::cout << "errorStage: " << errorStage << std::endl;

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
        else if (errorStage == 1) // this error stage stuff is not optimal, should pick the p guess from inital conditions
        {
            p = pow((sides[0]->a + sides[1]->a - 0.5 * (gammma - 1) * (sides[1]->u - sides[0]->u)) / (sides[0]->a / pow(sides[0]->p, (gammma - 1) / (2 * gammma)) + sides[1]->a / pow(sides[1]->p, (gammma - 1) / (2 * gammma))), (2 * gammma) / (gammma - 1));
//doesnt work because uR is too big
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
        else if (errorStage == 5)
        {
        }
        else if (errorStage == 6)
        {
            p = 1 * pow(10, -6);
        }
        else
        {
            std::cout << "rip bozo" << std::endl;
            return;
        }
        count = 0;
        while (iterate)
        {
            errno = 0;
            if (waveOutput)
                std::cout << "p* start: " << p << std::endl;

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
                if (waveOutput)
                    std::cout << "error occured \n"
                              << std::endl;

                errorStage++;
                break;
            }

            double f = fs[0] + fs[1] - sides[0]->u + sides[1]->u;
            double d_f = d_fs[0] + d_fs[1];

            change = f / d_f;

            p = p - change; // Update new estimate of p*

            if (waveOutput)
            {
                std::cout << "f: " << f << std::endl;
                std::cout << "d_f: " << d_f << std::endl;
                std::cout << "change: " << change << std::endl;
                std::cout << "p*: " << p << std::endl;
            }
            count++;

            if (waveOutput)
                std::cout << "End of iteration " << count << "\n"
                          << std::endl;

            // if (count == 5) {return 0;} // pause iteration and exit

            if (TOL >= 2 * fabs(change / (change + 2 * p))) // iteration limit (slightly different to notes as abs of entire rhs)
            {
                if (waveOutput)
                    std::cout << "iterate false" << std::endl;
                iterate = false;
            }
        }
    }
    if (domainOutput)
        std::cout << "loop ended count = " << count << ", reached errorStage " << errorStage << std::endl;

    u = 0.5 * (sides[0]->u + sides[1]->u) + 0.5 * (fs[1] - fs[0]); // u*

    if (u >= 0) // pick rho value depending on the side of the discontinuity
    {
        if (p > sides[0]->p)
        {
            rho = sides[0]->rho * (((p / sides[0]->p) + ((gammma - 1) / (gammma + 1))) / (((gammma - 1) / (gammma + 1)) * (p / sides[0]->p) + 1));
        }
        else
        {
            rho = sides[0]->rho * pow((p / sides[0]->p), 1 / gammma);
        }
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
    }

    aCalc();
    return;
}
