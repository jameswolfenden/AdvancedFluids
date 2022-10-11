#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>

const double gamma = 1.4;

class Point
{
public:
    double p, rho, u, a;

    Point(double p)
    {
        this->p = p;
    }

    Point(double p, double rho, double u)
    {
        this->p = p;
        this->rho = rho;
        this->u = u;
        calca();
    }

    // find speed of sound at point kinda pointless since * region has 2 different a's
    double calca()
    {
        a = sqrt((gamma * p) / rho);
        return a;
    }
};

int main()
{

    Point sides[] = {Point(1.0, 1.0, 0.0), Point(0.1, 0.125, 0)}; // Test case 1
    //Point sides[] = {Point(0.4, 1.0, -2.0), Point(0.4, 1.0, 2.0)}; // Test case 2

    Point star(0.5 * (sides[0].p + sides[1].p)); // Point at the star region

    double change;
    int count = 0;
    double fs[2];
    do 
    {

        double d_fs[2];

        for (int side = 0; side < 2; side++)
        {

            double A = 2 / ((gamma + 1) * sides[side].rho);
            double B = sides[side].p * (gamma - 1) / (gamma + 1);

            if (star.p > sides[side].p) // shock
            {
                fs[side] = (star.p - sides[side].p) * abs(pow(std::complex<double>(A / (star.p + B), 0), 0.5));
                d_fs[side] = abs(pow(std::complex<double>(A / (star.p + B), 0), 0.5)) * (1 - (star.p - sides[side].p) / (2 * (B + star.p)));
            }
            else // expansion
            {
                fs[side] = 2 * sides[side].a / (gamma - 1) * (abs(pow(std::complex<double>(star.p / sides[side].p, 0), (gamma - 1) / (2 * gamma))) - 1);
                d_fs[side] = 1 / (sides[side].p * sides[side].a) * abs(pow(std::complex<double>(star.p / sides[side].p, 0), -(gamma - 1) / (2 * gamma)));
            }
        }

        double f = fs[0] + fs[1] - sides[0].u + sides[1].u;
        double d_f = d_fs[0] + d_fs[1];

        change = f / d_f;

        star.p = star.p - change; // Update new estimate of p*

        std::cout << star.p << std::endl;

        count++;

        // if (count == 5) {return 0;}

    } while (!(pow(10, -6) >= 2 * fabs(change / (change + 2 * star.p)))); // iteration limit (slightly different to notes as abs of entire rhs)
    std::cout << "loop ended count = " << count << std::endl;
    std::cout << star.p << std::endl;
    std::cout << count << std::endl;

    star.u = 0.5 * (sides[0].u + sides[1].u) + 0.5 * (fs[1] - fs[0]); // u*

    std::cout << star.u << std::endl;

    double starRhos[2];
    double uShock, uHead, uTail;
    double rhoFan[100], uFan[100], pFan[100], u[100];
    for (int side = 0; side < 2; side++)
    {
        if (star.p > sides[side].p) // shock
        {
            // Find the rho* for both sides of the discontinuity for outside wave shock
            starRhos[side] = sides[side].rho * (((star.p / sides[side].p) + ((gamma - 1) / (gamma + 1))) / (((gamma - 1) / (gamma + 1)) * (star.p / sides[side].p) + 1));
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
            uShock = sides[side].u + S * sides[side].a * pow((gamma + 1) / (2 * gamma) * star.p / sides[side].p + (gamma - 1) / (2 * gamma), 0.5);
        }
        else // expansion
        {
            // rho*s in case outside wave is expansion
            starRhos[side] = sides[side].rho * pow((star.p / sides[side].p), 1 / gamma);
            int S;
            if (side == 0)
            {
                S = -1;
            }
            else
            {
                S = 1;
            }
            uHead = sides[side].u + S * sides[side].a; // head of expansion fan speed
            double aStarSide = sqrt((gamma * star.p) / starRhos[side]); // get speed of sound for the side but close to *
            uTail = star.u + S * aStarSide; // tail of expansion fan speed

            // Get 100 points for ploting properties within expansion fan
            for (int i = 0; i < 100; i++)
            {
                u[i] = (uHead - uTail) / 100 * i + uTail; // linear distribution of velocities within the fan
                rhoFan[i] = sides[side].rho * pow(2 / (gamma + 1) - S * (gamma - 1) / (gamma + 1) * (sides[side].u - u[i]) / sides[side].a, 2 / (gamma - 1));
                uFan[i] = 2 / (gamma + 1) * (-S * sides[side].a + ((gamma - 1) / 2) * sides[side].u + u[i]);
                pFan[i] = sides[side].p * pow(2 / (gamma + 1) - S * (gamma - 1) / (gamma + 1) * (sides[side].u - u[i]) / sides[side].a, (2 * gamma) / (gamma - 1));
            }
        }
    }

    // write vales to a csv for reading
    std::ofstream csvFile;
    csvFile.open("save2.csv");
    csvFile << "p*," << star.p << "\n";
    csvFile << "u*," << star.u << "\n";
    csvFile << "rho*L," << starRhos[0] << "\n";
    csvFile << "rho*R," << starRhos[1] << "\n";
    csvFile << "uShock," << uShock << "\n";
    csvFile << "uHead," << uHead << "\n";
    csvFile << "uTail," << uTail << "\n";
    csvFile << "\n";
    csvFile << "Fan stuff\n";
    csvFile << "u,rhoFan,uFan,pFan\n";
    for (int i = 0; i < 100; i++)
    {
        csvFile << u[i] << "," << rhoFan[i] << "," << uFan[i] << "," << pFan[i] << "\n";
    }
    csvFile.close();

    // write non-iterative values to csv for opening in matlab
    std::ofstream csvFileGraph1;
    csvFileGraph1.open("graph1.csv");
    csvFileGraph1 << star.p << "\n";
    csvFileGraph1 << star.u << "\n";
    csvFileGraph1 << starRhos[0] << "\n";
    csvFileGraph1 << starRhos[1] << "\n";
    csvFileGraph1 << uShock << "\n";
    csvFileGraph1 << uHead << "\n";
    csvFileGraph1 << uTail << "\n";
    csvFileGraph1.close();

    // write the fan values to csv for matlab
    std::ofstream csvFileGraph2;
    csvFileGraph2.open("graph2.csv");
    for (int i = 0; i < 100; i++)
    {
        csvFileGraph2 << u[i] << "," << rhoFan[i] << "," << uFan[i] << "," << pFan[i] << "\n";
    }
    csvFileGraph2.close();
}

// double fShock = (star.p - sides[side].p) * pow(A / (star.p + B), 0.5);
// double fExpansion = 2 * sides[side].a / (gamma - 1) * (pow(star.p / sides[side].p, (gamma - 1) / (2 * gamma)) - 1);

// double d_fShock = pow(A / (star.p + B), 0.5) * (1 - (star.p - sides[side].p) / (2 * (B + star.p)));
// double d_fExpansion = 1 / (sides[side].p * sides[side].a) * pow(star.p / sides[side].p, -(gamma - 1) / (2 * gamma));
