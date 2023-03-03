#include "Domain3D.h"
#include <iostream>
#include <fstream>
#include <algorithm>

void writeToDomainStreams(Domain3D *domain, std::vector<std::ofstream> &streams, int z)
{
    // write to all file streams in the vector
    for (int y = 0; y < (domain->yCellCount); y++) // write final cell pressures
    {
        for (int x = 0; x < domain->xCellCount; x++)
        {
            if (x == domain->xCellCount - 1)
            {
                streams[0] << domain->boxes[x][y][z].p;
                streams[1] << domain->boxes[x][y][z].rho;
                streams[2] << domain->boxes[x][y][z].u;
                streams[3] << domain->boxes[x][y][z].v;
                streams[4] << domain->boxes[x][y][z].w;
            }
            else
            {
                streams[0] << domain->boxes[x][y][z].p << ",";
                streams[1] << domain->boxes[x][y][z].rho << ",";
                streams[2] << domain->boxes[x][y][z].u << ",";
                streams[3] << domain->boxes[x][y][z].v << ",";
                streams[4] << domain->boxes[x][y][z].w << ",";
            }
        }
        streams[0] << "\n";
        streams[1] << "\n";
        streams[2] << "\n";
        streams[3] << "\n";
        streams[4] << "\n";
    }
}
void saveWallPressures(std::vector<Domain3D *> domains, std::ofstream &stream)
{
    for (int y = 1; y < domains[0]->yFaceCount; y++)
    {
        for (int i = 0; i < domains.size(); i++)
        {
            for (int z = 1; z < domains[i]->zFaceCount; z++)
            {
                if (i == domains.size() - 1 && z == domains[i]->zFaceCount - 1)
                {
                    stream << domains[i]->boxes[0][y][z].p;
                }
                else
                {
                    stream << domains[i]->boxes[0][y][z].p << ",";
                }
            }
        }
        stream << "\n";
    }
}

int main()
{
    int iterations = 50;
    double elapsedTime = 0;
    int domainCount = 10;
    double fridgeHeight = 2.0; // full size
    // double fridgeHeight = 0.20; // reduced size
    double pipeHeight = 0.20;
    double pipeWidth = 0.03;
    double leftFridgeWidth = 0.05;
    double rightFridgeWidth = 1.0 - leftFridgeWidth - pipeWidth; // maybe right? full size
    // double rightFridgeWidth = 0.12; // reduced size
    double pipeDepth = 0.03;
    // double frontFridgeDepth = 0.40; // reduced size
    // double backFridgeDepth = 0.40; // reduced size
    double frontFridgeDepth = (1 - pipeDepth) / 2; // full size
    double backFridgeDepth = (1 - pipeDepth) / 2;  // full size
    int xCellsPerMetre = 100;
    int yCellsPerMetre = 100;
    int zCellsPerMetre = 100;
    int pipeWidthCells = (int)(pipeWidth * xCellsPerMetre) + 2;
    int pipeHeightCells = (int)(pipeHeight * yCellsPerMetre) + 2;
    int leftFridgeWidthCells = (int)(leftFridgeWidth * xCellsPerMetre) + 2;
    int rightFridgeWidthCells = (int)(rightFridgeWidth * xCellsPerMetre) + 2;
    int fridgeHeightCells = (int)(fridgeHeight * yCellsPerMetre) + 2;
    int pipeDepthCells = (int)(pipeDepth * zCellsPerMetre) + 2;
    int frontFridgeDepthCells = (int)(frontFridgeDepth * zCellsPerMetre) + 2;
    int backFridgeDepthCells = (int)(backFridgeDepth * zCellsPerMetre) + 2;
    std::cout << "pipeWidthCells: " << pipeWidthCells << std::endl;
    std::cout << "pipeHeightCells: " << pipeHeightCells << std::endl;
    std::cout << "leftFridgeWidthCells: " << leftFridgeWidthCells << std::endl;
    std::cout << "rightFridgeWidthCells: " << rightFridgeWidthCells << std::endl;
    std::cout << "fridgeHeightCells: " << fridgeHeightCells << std::endl;
    std::cout << "pipeDepthCells: " << pipeDepthCells << std::endl;
    std::cout << "frontFridgeDepthCells: " << frontFridgeDepthCells << std::endl;
    std::cout << "backFridgeDepthCells: " << backFridgeDepthCells << std::endl;
    // important to take into account the ghost cells - the domain is 1 cell larger in each direction
    fridgeHeight = (double)fridgeHeightCells / yCellsPerMetre;
    pipeHeight = (double)pipeHeightCells / yCellsPerMetre;
    pipeWidth = (double)pipeWidthCells / xCellsPerMetre;
    leftFridgeWidth = (double)leftFridgeWidthCells / xCellsPerMetre;
    rightFridgeWidth = (double)rightFridgeWidthCells / xCellsPerMetre;
    pipeDepth = (double)pipeDepthCells / zCellsPerMetre;
    frontFridgeDepth = (double)frontFridgeDepthCells / zCellsPerMetre;
    backFridgeDepth = (double)backFridgeDepthCells / zCellsPerMetre;
    std::cout << "fridgeHeight: " << fridgeHeight << std::endl;
    std::cout << "pipeHeight: " << pipeHeight << std::endl;
    std::cout << "pipeWidth: " << pipeWidth << std::endl;
    std::cout << "leftFridgeWidth: " << leftFridgeWidth << std::endl;
    std::cout << "rightFridgeWidth: " << rightFridgeWidth << std::endl;
    std::cout << "pipeDepth: " << pipeDepth << std::endl;
    std::cout << "frontFridgeDepth: " << frontFridgeDepth << std::endl;
    std::cout << "backFridgeDepth: " << backFridgeDepth << std::endl;

    int z = (pipeDepthCells - 1) / 2; // centre of pipe
    std::cout << "z: " << z << std::endl;

    double fridgePressure = 1.01325;
    double pipePressure = 1.06;
    double fridgeDensity = 1.268;
    double pipeDensity = 1.37;

    std::vector<std::vector<bool>> domainsGhosts = {
        // order of sides is top of domain, bottom, left, right, front, back
        {true, false, false, false, false, false}, // centre fridge interior
        {false, true, true, true, true, true},     // pipe below fridge
        {true, true, true, false, false, false},   // left fridge interior (small)
        {true, true, false, true, false, false},   // right fridge interior (large)
        {true, true, false, false, true, false},   // centre fridge front
        {true, true, false, false, false, true},   // centre fridge back
        {true, true, true, false, true, false},    // left fridge front
        {true, true, true, false, false, true},    // left fridge back
        {true, true, false, true, true, false},    // right fridge front
        {true, true, false, true, false, true}     // right fridge back
    };
    std::vector<std::vector<Domain3D *>> domainsTransmissive;
    domainsTransmissive.resize(domainCount, std::vector<Domain3D *>(6));
    std::vector<Domain3D> domains = {
        Domain3D(pipeWidth, fridgeHeight, pipeDepth, pipeWidthCells, fridgeHeightCells, pipeDepthCells, domainsGhosts[0], fridgePressure, fridgeDensity, 0.0, 0.0, 0.0),                             // centre fridge interior
        Domain3D(pipeWidth, pipeHeight, pipeDepth, pipeWidthCells, pipeHeightCells, pipeDepthCells, domainsGhosts[1], pipePressure, pipeDensity, 0.0, 0.0, 0.0),                                     // pipe below fridge
        Domain3D(leftFridgeWidth, fridgeHeight, pipeDepth, leftFridgeWidthCells, fridgeHeightCells, pipeDepthCells, domainsGhosts[2], fridgePressure, fridgeDensity, 0.0, 0.0, 0.0),                 // left fridge interior (small)
        Domain3D(rightFridgeWidth, fridgeHeight, pipeDepth, rightFridgeWidthCells, fridgeHeightCells, pipeDepthCells, domainsGhosts[3], fridgePressure, fridgeDensity, 0.0, 0.0, 0.0),               // right fridge interior (large)
        Domain3D(pipeWidth, fridgeHeight, frontFridgeDepth, pipeWidthCells, fridgeHeightCells, frontFridgeDepthCells, domainsGhosts[4], fridgePressure, fridgeDensity, 0.0, 0.0, 0.0),               // centre fridge front
        Domain3D(pipeWidth, fridgeHeight, backFridgeDepth, pipeWidthCells, fridgeHeightCells, backFridgeDepthCells, domainsGhosts[5], fridgePressure, fridgeDensity, 0.0, 0.0, 0.0),                 // centre fridge back
        Domain3D(leftFridgeWidth, fridgeHeight, frontFridgeDepth, leftFridgeWidthCells, fridgeHeightCells, frontFridgeDepthCells, domainsGhosts[6], fridgePressure, fridgeDensity, 0.0, 0.0, 0.0),   // left fridge front
        Domain3D(leftFridgeWidth, fridgeHeight, backFridgeDepth, leftFridgeWidthCells, fridgeHeightCells, backFridgeDepthCells, domainsGhosts[7], fridgePressure, fridgeDensity, 0.0, 0.0, 0.0),     // left fridge back
        Domain3D(rightFridgeWidth, fridgeHeight, frontFridgeDepth, rightFridgeWidthCells, fridgeHeightCells, frontFridgeDepthCells, domainsGhosts[8], fridgePressure, fridgeDensity, 0.0, 0.0, 0.0), // right fridge front
        Domain3D(rightFridgeWidth, fridgeHeight, backFridgeDepth, rightFridgeWidthCells, fridgeHeightCells, backFridgeDepthCells, domainsGhosts[9], fridgePressure, fridgeDensity, 0.0, 0.0, 0.0)};  // right fridge back
    domainsTransmissive[0][1] = &domains[1];                                                                                                                                                         // bottom of centre fridge domain is transmissive to pipe domain
    domainsTransmissive[0][2] = &domains[2];                                                                                                                                                         // left of centre fridge domain is transmissive to left fridge domain
    domainsTransmissive[0][3] = &domains[3];                                                                                                                                                         // right of centre fridge domain is transmissive to right fridge domain
    domainsTransmissive[0][4] = &domains[4];                                                                                                                                                         // front of centre fridge domain is transmissive to front of centre fridge domain
    domainsTransmissive[0][5] = &domains[5];                                                                                                                                                         // back of centre fridge domain is transmissive to back of centre fridge domain
    domainsTransmissive[1][0] = &domains[0];                                                                                                                                                         // top of pipe domain is transmissive to centre fridge domain
    domainsTransmissive[2][3] = &domains[0];                                                                                                                                                         // right of left fridge domain is transmissive to centre fridge domain
    domainsTransmissive[2][4] = &domains[6];                                                                                                                                                         // front of left fridge domain is transmissive to front of left fridge domain
    domainsTransmissive[2][5] = &domains[7];                                                                                                                                                         // back of left fridge domain is transmissive to back of left fridge domain
    domainsTransmissive[3][2] = &domains[0];                                                                                                                                                         // left of right fridge domain is transmissive to centre fridge domain
    domainsTransmissive[3][4] = &domains[8];                                                                                                                                                         // front of right fridge domain is transmissive to front of right fridge domain
    domainsTransmissive[3][5] = &domains[9];                                                                                                                                                         // back of right fridge domain is transmissive to back of right fridge domain
    domainsTransmissive[4][5] = &domains[0];                                                                                                                                                         // back of front of centre fridge domain is transmissive to centre fridge domain
    domainsTransmissive[5][4] = &domains[0];                                                                                                                                                         // front of back of centre fridge domain is transmissive to centre fridge domain
    domainsTransmissive[6][5] = &domains[2];                                                                                                                                                         // back of front of left fridge domain is transmissive to left fridge domain
    domainsTransmissive[7][4] = &domains[2];                                                                                                                                                         // front of back of left fridge domain is transmissive to left fridge domain
    domainsTransmissive[8][5] = &domains[3];                                                                                                                                                         // back of front of right fridge domain is transmissive to right fridge domain
    domainsTransmissive[9][4] = &domains[3];                                                                                                                                                         // front of back of right fridge domain is transmissive to right fridge domain
    domainsTransmissive[4][2] = &domains[6];                                                                                                                                                         // left of front of centre fridge domain is transmissive to front of left fridge domain
    domainsTransmissive[4][3] = &domains[8];                                                                                                                                                         // right of front of centre fridge domain is transmissive to front of right fridge domain
    domainsTransmissive[5][2] = &domains[7];                                                                                                                                                         // left of back of centre fridge domain is transmissive to back of left fridge domain
    domainsTransmissive[5][3] = &domains[9];                                                                                                                                                         // right of back of centre fridge domain is transmissive to back of right fridge domain
    domainsTransmissive[6][3] = &domains[4];                                                                                                                                                         // right of front of left fridge domain is transmissive to front of centre fridge domain
    domainsTransmissive[7][3] = &domains[5];                                                                                                                                                         // right of back of left fridge domain is transmissive to back of centre fridge domain
    domainsTransmissive[8][2] = &domains[4];                                                                                                                                                         // left of front of right fridge domain is transmissive to front of centre fridge domain
    domainsTransmissive[9][2] = &domains[5];                                                                                                                                                         // left of back of right fridge domain is transmissive to back of centre fridge domain

    // 2d vector of file streams for each domain
    std::vector<std::vector<std::ofstream>> domainStreams;
    // create a vector of file streams for each domain and add to the 2d vector
    for (int domain = 0; domain < 4; domain++)
    {
        std::vector<std::ofstream> domainStream;
        for (int i = 0; i < 5; i++)
        {
            std::ofstream stream;
            domainStream.push_back(std::move(stream));
        }
        domainStreams.push_back(std::move(domainStream));
    }
    for (int domain = 0; domain < 4; domain++)
    {
        domainStreams[domain][0].open("./resultsdomain/domain" + std::to_string(domain) + "p.csv");
        domainStreams[domain][1].open("./resultsdomain/domain" + std::to_string(domain) + "rho.csv");
        domainStreams[domain][2].open("./resultsdomain/domain" + std::to_string(domain) + "u.csv");
        domainStreams[domain][3].open("./resultsdomain/domain" + std::to_string(domain) + "v.csv");
        domainStreams[domain][4].open("./resultsdomain/domain" + std::to_string(domain) + "w.csv");
    }

    // save initial state
    for (int domain = 0; domain < 4; domain++)
    {
        writeToDomainStreams(&domains[domain], domainStreams[domain], z);
    }

    std::ofstream wallPressuresStream;
    wallPressuresStream.open("./resultsdomain/wallPressures.csv");
    std::vector<Domain3D *> wallDomains = {&domains[6], &domains[2], &domains[7]};
    saveWallPressures(wallDomains, wallPressuresStream);

    bool run = true;
    for (int i = 1; i < iterations + 1; i++) // iterate updating the domain's cells
    {
        if (run)
        {
            std::cout << "Starting iteration " << i << std::endl;
            double timeStep = std::min({domains[0].timeStep(), domains[1].timeStep(), domains[2].timeStep(), domains[3].timeStep()}); // find smallest allowable timestep from all domains
            for (int domain = 0; domain < domainCount; domain++)
            {
                if (!domains[domain].updateBoxes(domainsTransmissive[domain], timeStep)) // update each domain by the allowable timestep, pointer to the transmissive domains passed
                {
                    std::cout << "Error updating domain " << domain << std::endl;
                    run = false;
                    break;
                }
            }
            elapsedTime += (timeStep * 2); // not entirely sure its x2 but i think so as it iterates 2 with xyyx

            // save state
            for (int domain = 0; domain < 4; domain++)
            {
                writeToDomainStreams(&domains[domain], domainStreams[domain], z);
            }
            saveWallPressures(wallDomains, wallPressuresStream);
        }
        else
        {
            std::ofstream errorStream;
            errorStream.open("./resultsdomain/error.txt");
            errorStream << "Error updating domain";
            errorStream.close();
            break;
        }
    }
    // close file streams
    for (int domain = 0; domain < 4; domain++)
    {
        for (int stream = 0; stream < 4; stream++)
        {
            domainStreams[domain][stream].close();
        }
    }
    wallPressuresStream.close();
    // save number of iterations to file
    std::ofstream iterationsStream;
    iterationsStream.open("./resultsdomain/iterations.csv");
    iterationsStream << iterations;
    iterationsStream.close();

    std::cout << elapsedTime << std::endl; // output elapsed time (of simulation time) to console
}
