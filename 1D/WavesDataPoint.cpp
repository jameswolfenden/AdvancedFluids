#include "WavesDataPoint.h"
#include <cmath>
#include "../fluidConsts.h"
#include <iostream>

void WavesDataPoint::waveData(Point *sides[])
{
    for (int side = 0; side < 2; side++)
    {
        if (p > sides[side]->p) // shock
        {
            // Find the rho* for both sides of the discontinuity for outside wave shock
            rhos[side] = sides[side]->rho * (((p / sides[side]->p) + ((gammma - 1) / (gammma + 1))) / (((gammma - 1) / (gammma + 1)) * (p / sides[side]->p) + 1));
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
            uShock = sides[side]->u + S * sides[side]->a * pow((gammma + 1) / (2 * gammma) * p / sides[side]->p + (gammma - 1) / (2 * gammma), 0.5);
        }
        else // expansion
        {
            // rho*s in case outside wave is expansion
            rhos[side] = sides[side]->rho * pow((p / sides[side]->p), 1 / gammma);
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
            double aStarSide = sqrt((gammma * p) / rhos[side]); // get speed of sound for the side but close to *
            uTail = u + S * aStarSide;                         // tail of expansion fan speed

            // Get 100 points for ploting properties within expansion fan
            for (int i = 0; i < 100; i++)
            {
                us[i] = (uHead - uTail) / 100 * i + uTail; // linear distribution of velocities within the fan
                rhoFan[i] = sides[side]->rho * pow(2 / (gammma + 1) - S * (gammma - 1) / (gammma + 1) * (sides[side]->u - us[i]) / sides[side]->a, 2 / (gammma - 1));
                uFan[i] = 2 / (gammma + 1) * (-S * sides[side]->a + ((gammma - 1) / 2) * sides[side]->u + us[i]);
                pFan[i] = sides[side]->p * pow(2 / (gammma + 1) - S * (gammma - 1) / (gammma + 1) * (sides[side]->u - us[i]) / sides[side]->a, (2 * gammma) / (gammma - 1));
            }
        }
    }
    return;
}
