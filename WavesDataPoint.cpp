#include "WavesDataPoint.h"
#include <cmath>
#include "fluidConsts.h"
#include <iostream>

void WavesDataPoint::findStar(Point *sides[])
{

    double change;
    int count = 0;
    int errorStage = 0;
    bool iterate = true;
    double fs[2];
    double d_fs[2];

    while (iterate)
    {
        std::cout << "errorStage: " << errorStage << std::endl;

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
        else if (errorStage == 3)
        {
        }
        else if (errorStage == 4)
        {
        }
        else if (errorStage == 5)
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
            std::cout << "p* start: " << p << std::endl;

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
                std::cout << "error occured \n"
                          << std::endl;

                errorStage++;
                break;
            }

            double f = fs[0] + fs[1] - sides[0]->u + sides[1]->u;
            double d_f = d_fs[0] + d_fs[1];

            change = f / d_f;

            p = p - change; // Update new estimate of p*

            std::cout << "f: " << f << std::endl;
            std::cout << "d_f: " << d_f << std::endl;
            std::cout << "change: " << change << std::endl;
            std::cout << "p*: " << p << std::endl;

            count++;

            std::cout << "End of iteration " << count << "\n"
                      << std::endl;

            // if (count == 5) {return 0;} // pause iteration and exit

            if (TOL >= 2 * fabs(change / (change + 2 * p))) // iteration limit (slightly different to notes as abs of entire rhs)
            {
                std::cout << "iterate false" << std::endl;
                iterate = false;
            }
        }
    }
    std::cout << "loop ended count = " << count << ", reached errorStage " << errorStage << std::endl;

    u = 0.5 * (sides[0]->u + sides[1]->u) + 0.5 * (fs[1] - fs[0]); // u*
    return;
}

void WavesDataPoint::waveData(Point *sides[])
{
    for (int side = 0; side < 2; side++)
    {
        if (p > sides[side]->p) // shock
        {
            // Find the rho* for both sides of the discontinuity for outside wave shock
            rhos[side] = sides[side]->rho * (((p / sides[side]->p) + ((gamma - 1) / (gamma + 1))) / (((gamma - 1) / (gamma + 1)) * (p / sides[side]->p) + 1));
            int S;
            if (side == 0)
            {
                S = -1;
            }
            else
            {
                S = 1;
            }
            // Find the shock speed
            uShock = sides[side]->u + S * sides[side]->a * pow((gamma + 1) / (2 * gamma) * p / sides[side]->p + (gamma - 1) / (2 * gamma), 0.5);
        }
        else // expansion
        {
            // rho*s in case outside wave is expansion
            rhos[side] = sides[side]->rho * pow((p / sides[side]->p), 1 / gamma);
            int S;
            if (side == 0)
            {
                S = -1;
            }
            else
            {
                S = 1;
            }
            uHead = sides[side]->u + S * sides[side]->a;       // head of expansion fan speed
            double aStarSide = sqrt((gamma * p) / rhos[side]); // get speed of sound for the side but close to *
            uTail = u + S * aStarSide;                         // tail of expansion fan speed

            // Get 100 points for ploting properties within expansion fan
            for (int i = 0; i < 100; i++)
            {
                us[i] = (uHead - uTail) / 100 * i + uTail; // linear distribution of velocities within the fan
                rhoFan[i] = sides[side]->rho * pow(2 / (gamma + 1) - S * (gamma - 1) / (gamma + 1) * (sides[side]->u - us[i]) / sides[side]->a, 2 / (gamma - 1));
                uFan[i] = 2 / (gamma + 1) * (-S * sides[side]->a + ((gamma - 1) / 2) * sides[side]->u + us[i]);
                pFan[i] = sides[side]->p * pow(2 / (gamma + 1) - S * (gamma - 1) / (gamma + 1) * (sides[side]->u - us[i]) / sides[side]->a, (2 * gamma) / (gamma - 1));
            }
        }
    }
    return;
}
