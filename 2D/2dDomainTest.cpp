#include "Domain2D.h"
#include <iostream>
#include <fstream>
#include <algorithm>

void saveDomain(Domain2D *domain, int num, int iteration) // each domain is saved seperatly 
{
    std::ofstream pCells, rhoCells, uCells, vCells;
    pCells.open("./results/"+std::to_string(num)+"pCells" + std::to_string(iteration) + ".csv");
    rhoCells.open("./results/"+std::to_string(num)+"rhoCells" + std::to_string(iteration) + ".csv");
    uCells.open("./results/"+std::to_string(num)+"uCells" + std::to_string(iteration) + ".csv");
    vCells.open("./results/"+std::to_string(num)+"vCells" + std::to_string(iteration) + ".csv");
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

    std::vector<std::vector<bool>> domainsGhosts = { // order of sides is top of domain, bottom, left, right
        {false, true, false, false}, // centre bottom
        {true, false, true, true}, // top
        {true, true, true, false}, // left
        {true, true, false, true}}; // right
    std::vector<std::vector<Domain2D *>> domainsTransmissive;
    domainsTransmissive.resize(domainCount, std::vector<Domain2D *>(4));
    std::vector<Domain2D> domains = {
        Domain2D(1, 1, xCellCount, yCellCount, domainsGhosts[0], 0.1, 0.125, 0.0, 0.0), // centre bottom
        Domain2D(1, 1, xCellCount, yCellCount, domainsGhosts[1], 1, 1, 0.0, 0.0), // top
        Domain2D(1, 1, xCellCount, yCellCount, domainsGhosts[2], 0.1, 0.125, 0.0, 0.0), // left
        Domain2D(1, 1, xCellCount, yCellCount, domainsGhosts[3], 0.1, 0.125, 0.0, 0.0)}; //right
    domainsTransmissive[0][0] = &domains[1]; // top of centre bottom domain is transmissive to top domain
    domainsTransmissive[0][2] = &domains[2]; // left of centre bottom domain is transmissive to left domain
    domainsTransmissive[0][3] = &domains[3]; // right of centre bottom domain is transmissive to right domain
    domainsTransmissive[1][1] = &domains[0]; // bottom of top domain is transmissive to bottom domain
    domainsTransmissive[2][3] = &domains[0]; // right of left domain is transmissive to centre bottom domain
    domainsTransmissive[3][2] = &domains[0]; // left of right domain is transmissive to centre bottom domain

    for (int i = 0; i < iterations; i++) // iterate updating the domain's cells
    {
        std::cout << "Starting iteration " << i << std::endl;
        double timeStep = std::min({domains[0].timeStep(), domains[1].timeStep(), domains[2].timeStep(), domains[3].timeStep()}); // find smallest allowable timestep from all domains
        for (int domain = 0; domain < domainCount; domain++)
        {
            domains[domain].updateCells(domainsTransmissive[domain], timeStep); // update each domain by the allowable timestep, pointer to the transmissive domains passed
        }
        elapsedTime += (timeStep * 2); // not entirely sure its x2 but i think so as it iterates 2 with xyyx
    

    for (int domain = 0; domain < domainCount; domain++)
    {
        saveDomain(&domains[domain], domain, i);
    }
    }

    std::cout << elapsedTime << std::endl; // output elapsed time (of simulation time) to console
}
