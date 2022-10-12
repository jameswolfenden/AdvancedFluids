#include <iostream>
#include <cmath>
#include <fstream>
#include <errno.h>
#include "Point.h"
#include "fluidConsts.h"
#include "WavesDataPoint.h"

void saveToCSV(WavesDataPoint star, Point *sides[2])
{

    // write vales to a csv for reading
    std::ofstream csvFile;
    csvFile.open("save2.csv");
    csvFile << "p*," << star.p << "\n";
    csvFile << "u*," << star.u << "\n";
    csvFile << "rho*L," << star.rhos[0] << "\n";
    csvFile << "rho*R," << star.rhos[1] << "\n";
    csvFile << "uShock," << star.uShock << "\n";
    csvFile << "uHead," << star.uHead << "\n";
    csvFile << "uTail," << star.uTail << "\n";
    csvFile << "\n";
    csvFile << "Fan stuff\n";
    csvFile << "u,rhoFan,uFan,pFan\n";
    for (int i = 0; i < 100; i++)
    {
        csvFile << star.us[i] << "," << star.rhoFan[i] << "," << star.uFan[i] << "," << star.pFan[i] << "\n";
    }
    csvFile.close();

    // write non-iterative values to csv for opening in matlab
    std::ofstream csvFileGraph1;
    csvFileGraph1.open("graph1.csv");
    csvFileGraph1 << star.p << "\n";
    csvFileGraph1 << star.u << "\n";
    csvFileGraph1 << star.rhos[0] << "\n";
    csvFileGraph1 << star.rhos[1] << "\n";
    csvFileGraph1 << star.uShock << "\n";
    csvFileGraph1 << star.uHead << "\n";
    csvFileGraph1 << star.uTail << "\n";
    csvFileGraph1 << sides[0]->p << "\n";
    csvFileGraph1 << sides[0]->rho << "\n";
    csvFileGraph1 << sides[0]->u << "\n";
    csvFileGraph1 << sides[1]->p << "\n";
    csvFileGraph1 << sides[1]->rho << "\n";
    csvFileGraph1 << sides[1]->u << "\n";
    csvFileGraph1.close();

    // write the fan values to csv for matlab
    std::ofstream csvFileGraph2;
    csvFileGraph2.open("graph2.csv");
    for (int i = 0; i < 100; i++)
    {
        csvFileGraph2 << star.us[i] << "," << star.rhoFan[i] << "," << star.uFan[i] << "," << star.pFan[i] << "\n";
    }
    csvFileGraph2.close();
}

int main()
{

    Point left(1.0, 1.0, 0.0);
    Point right(0.1, 0.125, 0);

    // Point left(0.4, 1.0, -2.0);
    // Point right(0.4, 1.0, 2.0);

    Point *sides[2];
    sides[0] = &left;
    sides[1] = &right;

    WavesDataPoint star;
    star.findStar(sides);
    star.waveData(sides);
    saveToCSV(star, sides);
}