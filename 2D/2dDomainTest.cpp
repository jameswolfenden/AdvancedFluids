#include "Domain2D.h"
#include <iostream>
#include <fstream>
#include <algorithm>

void saveDomain(Domain2D *domain, int num, int iteration) // each domain is saved seperatly
{
    std::ofstream pCells, rhoCells, uCells, vCells;
    pCells.open("./results/" + std::to_string(num) + "pCells" + std::to_string(iteration) + ".csv");
    rhoCells.open("./results/" + std::to_string(num) + "rhoCells" + std::to_string(iteration) + ".csv");
    uCells.open("./results/" + std::to_string(num) + "uCells" + std::to_string(iteration) + ".csv");
    vCells.open("./results/" + std::to_string(num) + "vCells" + std::to_string(iteration) + ".csv");
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
    int iterations = 50;
    double elapsedTime = 0;
    int domainCount = 4;
   // double fridgeHeight = 2.02; // full size
    double fridgeHeight = 0.72; // reduced size
    double pipeHeight = 0.77;
    double pipeWidth = 0.05;
    double leftFridgeWidth = 0.07;
    //double rightFridgeWidth = 1.06 - leftFridgeWidth - pipeWidth; // maybe right? full size
    double rightFridgeWidth = 0.36 - leftFridgeWidth - pipeWidth; // reduced size
    int xCellsPerMetre = 200;
    int yCellsPerMetre = 200;
    int pipeWidthCells = (int) (pipeWidth*xCellsPerMetre);
    int pipeHeightCells = (int) (pipeHeight*yCellsPerMetre);
    int leftFridgeWidthCells = (int) (leftFridgeWidth*xCellsPerMetre);
    int rightFridgeWidthCells = (int) (rightFridgeWidth*xCellsPerMetre);
    int fridgeHeightCells = (int) (fridgeHeight*yCellsPerMetre);
    std::cout << "pipeWidthCells: " << pipeWidthCells << std::endl;
    std::cout << "pipeHeightCells: " << pipeHeightCells << std::endl;
    std::cout << "leftFridgeWidthCells: " << leftFridgeWidthCells << std::endl;
    std::cout << "rightFridgeWidthCells: " << rightFridgeWidthCells << std::endl;
    std::cout << "fridgeHeightCells: " << fridgeHeightCells << std::endl;
    // important to take into account the ghost cells - the domain is 1 cell larger in each direction



    std::vector<std::vector<bool>> domainsGhosts = {                             // order of sides is top of domain, bottom, left, right
                                                    {true, false, false, false}, // centre fridge interior
                                                    {false, true, true, true},   // pipe below fridge
                                                    {true, true, true, false},   // left fridge interior (small)
                                                    {true, true, false, true}};  // right fridge interior (large)
    std::vector<std::vector<Domain2D *>> domainsTransmissive;
    domainsTransmissive.resize(domainCount, std::vector<Domain2D *>(4));
    std::vector<Domain2D> domains = {
        Domain2D(pipeWidth, fridgeHeight, pipeWidthCells, fridgeHeightCells, domainsGhosts[0], 1, 1.3, 0.0, 0.0),  // centre fridge interior
        Domain2D(pipeWidth, pipeHeight, pipeWidthCells, pipeHeightCells, domainsGhosts[1], 1.1, 1.45, 0.0, 0.0),        // pipe below fridge
        Domain2D(leftFridgeWidth, fridgeHeight, leftFridgeWidthCells, fridgeHeightCells, domainsGhosts[2], 1, 1.3, 0.0, 0.0),  // left fridge interior (small)
        Domain2D(rightFridgeWidth, fridgeHeight, rightFridgeWidthCells, fridgeHeightCells, domainsGhosts[3], 1, 1.3, 0.0, 0.0)}; // right fridge interior (large)
    domainsTransmissive[0][1] = &domains[1];                                             // bottom of centre fridge domain is transmissive to pipe domain
    domainsTransmissive[0][2] = &domains[2];                                             // left of centre fridge domain is transmissive to left fridge domain
    domainsTransmissive[0][3] = &domains[3];                                             // right of centre fridge domain is transmissive to right fridge domain
    domainsTransmissive[1][0] = &domains[0];                                             // top of pipe domain is transmissive to centre fridge domain
    domainsTransmissive[2][3] = &domains[0];                                             // right of left fridge domain is transmissive to centre fridge domain
    domainsTransmissive[3][2] = &domains[0];                                             // left of right fridge domain is transmissive to centre fridge domain

    for (int domain = 0; domain < domainCount; domain++)
    {
        saveDomain(&domains[domain], domain, 0);
    }

    for (int i = 1; i < iterations+1; i++) // iterate updating the domain's cells
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
