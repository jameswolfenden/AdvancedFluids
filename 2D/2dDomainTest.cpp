#include "Domain2D.h"
#include <iostream>
#include <fstream>
#include <algorithm>

void saveDomain(Domain2D *domain, int num)
{
    std::ofstream pCells, rhoCells, uCells, vCells;
    pCells.open("pCells" + std::to_string(num) + ".csv");
    rhoCells.open("rhoCells" + std::to_string(num) + ".csv");
    uCells.open("uCells" + std::to_string(num) + ".csv");
    vCells.open("vCells" + std::to_string(num) + ".csv");
    for (int y = 0; y < (domain->yCellCount); y++) // write final cell pressures
    {
        for (int x = 0; x < domain->xCellCount; x++)
        {
            pCells << domain->cells[x][y].p << ",";
            rhoCells << domain->cells[x][y].rho << ",";
            uCells << domain->cells[x][y].u << ",";
            vCells << domain->cells[x][y].v << ",";
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
    int xCellCount = 100;
    int yCellCount = 100;
    int iterations = 100;
    double elapsedTime = 0;
    int domainCount = 4;

    std::vector<std::vector<bool>> domainsGhosts = {
        {false, true, false, false},
        {true, false, true, true},
        {true, true, true, false},
        {true, true, false, true}};
    std::vector<std::vector<Domain2D *>> domainsTransmissive;
    domainsTransmissive.resize(domainCount, std::vector<Domain2D *>(4));
    std::vector<Domain2D> domains = {
        Domain2D(1, 1, xCellCount, yCellCount, domainsGhosts[0], 0.1, 0.125, 0.0, 0.0),
        Domain2D(1, 1, xCellCount, yCellCount, domainsGhosts[1], 1, 1, 0.0, 0.0),
        Domain2D(1, 1, xCellCount, yCellCount, domainsGhosts[2], 0.1, 0.125, 0.0, 0.0),
        Domain2D(1, 1, xCellCount, yCellCount, domainsGhosts[3], 0.1, 0.125, 0.0, 0.0)};
    domainsTransmissive[0][0] = &domains[1];
    domainsTransmissive[0][2] = &domains[2];
    domainsTransmissive[0][3] = &domains[3];
    domainsTransmissive[1][1] = &domains[0];
    domainsTransmissive[2][3] = &domains[0];
    domainsTransmissive[3][2] = &domains[0];

    for (int i = 0; i < iterations; i++) // iterate updating the domain's cells
    {
        std::cout << "Starting iteration " << i << std::endl;
        double timeStep = std::min({domains[0].timeStep(), domains[1].timeStep(), domains[2].timeStep(), domains[3].timeStep()});
        for (int domain = 0; domain < domainCount; domain++)
        {
            domains[domain].updateCells(domainsTransmissive[domain], timeStep);
        }
        elapsedTime += (timeStep * 2); // not entirely sure its x2 but i think so as it iterates 2 with xyyx
    }

    for (int domain = 0; domain < domainCount; domain++)
    {
        saveDomain(&domains[domain], domain);
    }

    std::cout << elapsedTime << std::endl; // output elapsed time (of simulation time) to console
}
