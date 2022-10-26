#include "Domain2D.h"
#include <iostream>
#include <fstream>

int main()
{
    int xCellCount = 64;
    int yCellCount = 64;
    int iterations = 10;
    double elapsedTime = 0;

    std::vector<bool> domain1Ghosts = {false, true, true, true};
    std::vector<bool> domain2Ghosts = {true, false, true, true};
    std::vector<Domain2D*> domain1Transmissive;
    std::vector<Domain2D*> domain2Transmissive;
    domain1Transmissive.resize(4);
    domain2Transmissive.resize(4);

    Domain2D domain1(1, 1, xCellCount, yCellCount, domain1Ghosts, 0.1, 0.125, 0.0, 0.0); // create 1x1 physical size domain
    Domain2D domain2(1, 1, xCellCount, yCellCount, domain2Ghosts, 1.0, 1.0, 0.0, 0.0); // create 1x1 physical size domain

    domain1Transmissive[0] = &domain2;
    domain2Transmissive[1] = &domain1;

    std::ofstream pStartCells, rhoStartCells, uStartCells, vStartCells, pEndCells, rhoEndCells, uEndCells, vEndCells; // create file writers for the output

    pStartCells.open("pStartCells.csv");
    rhoStartCells.open("rhoStartCells.csv");
    uStartCells.open("uStartCells.csv");
    vStartCells.open("vStartCells.csv");
    for (int y = 0; y < (yCellCount); y++) // write inital start cell pressures
    {
        for (int x = 0; x < xCellCount; x++)
        {
            pStartCells << domain1.cells[x][y].p << ",";
            rhoStartCells << domain1.cells[x][y].rho << ",";
            uStartCells << domain1.cells[x][y].u << ",";
            vStartCells << domain1.cells[x][y].v << ",";
        }
        pStartCells << "\n";
        rhoStartCells << "\n";
        uStartCells << "\n";
        vStartCells << "\n";
    }
    pStartCells.close();
    rhoStartCells.close();
    uStartCells.close();
    vStartCells.close();

    for (int i = 0; i < iterations; i++) // iterate updating the domain's cells
    {
        double timeStep;
        if (domain1.timeStep() < domain2.timeStep()) {
            timeStep = domain1.timeStep();
        }
        else
        {
            timeStep = domain2.timeStep();
        }
        domain1.updateCells(domain1Transmissive, timeStep);
        domain2.updateCells(domain2Transmissive, timeStep);
        elapsedTime += (timeStep*2); // not entirely sure its x2 but i think so as it iterates 2 with xyyx
    }

    pEndCells.open("pEndCells.csv");
    rhoEndCells.open("rhoEndCells.csv");
    uEndCells.open("uEndCells.csv");
    vEndCells.open("vEndCells.csv");
    for (int y = 0; y < (yCellCount); y++) // write final cell pressures
    {
        for (int x = 0; x < xCellCount; x++)
        {
            pEndCells << domain1.cells[x][y].p << ",";
            rhoEndCells << domain1.cells[x][y].rho << ",";
            uEndCells << domain1.cells[x][y].u << ",";
            vEndCells << domain1.cells[x][y].v << ",";
        }
        pEndCells << "\n";
        rhoEndCells << "\n";
        uEndCells << "\n";
        vEndCells << "\n";
    }
    pEndCells.close();
    rhoEndCells.close();
    uEndCells.close();
    vEndCells.close();

    std::cout << elapsedTime << std::endl; // output elapsed time (of simulation time) to console
}