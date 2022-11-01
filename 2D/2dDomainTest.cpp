#include "Domain2D.h"
#include <iostream>
#include <fstream>
#include <algorithm>


void saveDomain(Domain2D domain, int num) {
    std::ofstream pCells, rhoCells, uCells, vCells;
    pCells.open("pCells"+std::to_string(num)+".csv");
    rhoCells.open("rhoCells"+std::to_string(num)+".csv");
    uCells.open("uCells"+std::to_string(num)+".csv");
    vCells.open("vCells"+std::to_string(num)+".csv");
    for (int y = 0; y < (domain.yCellCount); y++) // write final cell pressures
    {
        for (int x = 0; x < domain.xCellCount; x++)
        {
            pCells << domain.cells[x][y].p << ",";
            rhoCells << domain.cells[x][y].rho << ",";
            uCells << domain.cells[x][y].u << ",";
            vCells << domain.cells[x][y].v << ",";
        }
        pCells << "\n";
        rhoCells << "\n";
        uCells << "\n";
        vCells << "\n";
    }
    pCells.close();
    rhoCells.close();
    uCells.close();
    vCells.close();
}

int main()
{
    int xCellCount = 64;
    int yCellCount = 64;
    int iterations = 40;
    double elapsedTime = 0;

    std::vector<bool> domain1Ghosts = {false, true, false, false};
    std::vector<bool> domain2Ghosts = {true, false, true, true};
    std::vector<bool> domain3Ghosts = {true, true, true, false};
    std::vector<bool> domain4Ghosts = {true, true, false, true};
    std::vector<Domain2D*> domain1Transmissive;
    std::vector<Domain2D*> domain2Transmissive;
    std::vector<Domain2D*> domain3Transmissive;
    std::vector<Domain2D*> domain4Transmissive;
    domain1Transmissive.resize(4);
    domain2Transmissive.resize(4);
    domain3Transmissive.resize(4);
    domain4Transmissive.resize(4);
    Domain2D domain1(1, 1, xCellCount, yCellCount, domain1Ghosts, 0.1, 0.125, 0.0, 0.0); // create 1x1 physical size domain
    Domain2D domain2(1, 1, xCellCount, yCellCount, domain2Ghosts, 1, 1, 0.0, 0.0); // create 1x1 physical size domain
    Domain2D domain3(1, 1, xCellCount, yCellCount, domain3Ghosts, 0.1, 0.125, 0.0, 0.0); // create 1x1 physical size domain
    Domain2D domain4(1, 1, xCellCount, yCellCount, domain4Ghosts, 0.1, 0.125, 0.0, 0.0); // create 1x1 physical size domain
    domain1Transmissive[0] = &domain2;
    domain1Transmissive[2] = &domain3;
    domain1Transmissive[3] = &domain4;
    domain2Transmissive[1] = &domain1;
    domain3Transmissive[3] = &domain1;
    domain4Transmissive[2] = &domain1;

    for (int i = 0; i < iterations; i++) // iterate updating the domain's cells
    {
        double timeStep = std::min({domain1.timeStep(), domain2.timeStep(),domain3.timeStep(), domain4.timeStep()});
        domain1.updateCells(domain1Transmissive, timeStep);
        domain2.updateCells(domain2Transmissive, timeStep);
        domain3.updateCells(domain3Transmissive, timeStep);
        domain4.updateCells(domain4Transmissive, timeStep);
        elapsedTime += (timeStep*2); // not entirely sure its x2 but i think so as it iterates 2 with xyyx
    }

    

    saveDomain(domain1,1);
    saveDomain(domain2,2);
    saveDomain(domain3,3);
    saveDomain(domain4,4);
    

    std::cout << elapsedTime << std::endl; // output elapsed time (of simulation time) to console
}
