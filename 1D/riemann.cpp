#include <iostream>
#include <cmath>
#include <fstream>
#include <errno.h>
#include "Point.h"
#include "../fluidConsts.h"
#include "WavesDataPoint.h"
#include "Domain.h"
#include <vector>
#include <chrono>

void saveToCSV(WavesDataPoint star, Point *sides[2]) // for saving the results of the riemann tests
{

    // write vales to a csv for reading
    std::ofstream csvFile;
    csvFile.open("save2.csv");
    csvFile << "p*," << star.p << "\n";
    csvFile << "u*," << star.u << "\n";
    csvFile << "rho*," << star.rho << "\n";
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

void sodTests(int caseTest) // run a sod test case and get the results to csv
{
    Point left;
    Point right;
    Point *_left = &left;
    Point *_right = &right;
    if (caseTest == 1)
    {
        _left->updatePrimatives(1.0, 1.0, 0.0);
        _right->updatePrimatives(0.1, 0.125, 0);
    }
    else if (caseTest == 2)
    {
        _left->updatePrimatives(0.4, 1.0, -2.0);
        _right->updatePrimatives(0.4, 1.0, 2.0);
    }
    else if (caseTest == 3)
    {
        _left->updatePrimatives(1000.0, 1.0, 0.0);
        _right->updatePrimatives(0.01, 1.0, 0.0);
    }
    else if (caseTest == 4)
    {
        _left->updatePrimatives(0.01, 1.0, 0.0);
        _right->updatePrimatives(100.0, 1.0, 0.0);
    }
    else if (caseTest == 5)
    {
        _left->updatePrimatives(460.894, 5.99924, 19.5975);
        _right->updatePrimatives(46.0950, 5.99242, -6.19633);
    }
    else
    {
        _left->updatePrimatives(0.76657, 1.24179, 0.829351);
        _right->updatePrimatives(0.000043548100000000, 0.000486327, -0.0234923);
    }

    Point *sides[2];
    sides[0] = _left;
    sides[1] = _right;

    WavesDataPoint star;
    star.findStar(sides);
    star.waveData(sides);
    saveToCSV(star, sides);
}

void domainTest() // run a test on a domain (1d)
{
    std::cout << "domaintest" << std::endl;
    int pointsSize = 1000;
    Domain domain(1, pointsSize); // create new domain of physical size 1
    int iterations = 7500;
    // create 2d vectors of rho,u,p as arrays don't can't be this large
    // Should change to vector of points or get rid of entirely
    std::vector<std::vector<double>> rhos;
    rhos.resize(iterations, std::vector<double>(pointsSize, 1)); // resize the vectors first to avoid resizing in loops
    std::vector<std::vector<double>> us;
    us.resize(iterations, std::vector<double>(pointsSize, 0));
    std::vector<std::vector<double>> ps;
    ps.resize(iterations, std::vector<double>(pointsSize, 0));
    std::vector<double> ts; // 1d vector of time
    ts.resize(iterations, 0);
    std::cout << "done creating arrays, starting iteration..." << std::endl;

    std::chrono::steady_clock::time_point beginIteration = std::chrono::steady_clock::now();
    std::ofstream rhoCSV, uCSV, pCSV, tCSV; // create 4 file writers for the output
    rhoCSV.open("rho.csv");
    uCSV.open("u.csv");
    pCSV.open("p.csv");
    tCSV.open("t.csv");
    for (int i = 0; i < iterations; i++)
    {
        domain.updatePoints();                                     // update the points in the domain
        for (int point = 0; point < domain.points.size(); point++) // add the new point properties to the vectors
        {
            // rhos[i][point] = domain.points[point].rho;
            // us[i][point] = domain.points[point].u;
            // ps[i][point] = domain.points[point].p;
            rhoCSV << domain.points[point].rho;
            uCSV << domain.points[point].u;
            pCSV << domain.points[point].p;
            if (point + 1 < domain.points.size())
            {
                rhoCSV << ",";
                uCSV << ",";
                pCSV << ",";
            }
        }
        // ts[i] = domain.elapsedTime; // put the elapsed time into the vector for the current iteration
        rhoCSV << "\n"; // use n not endl as much faster
        uCSV << "\n";
        pCSV << "\n";
        tCSV << domain.elapsedTime << "\n";
    }
    rhoCSV.close();
    uCSV.close();
    pCSV.close();
    tCSV.close();
    std::chrono::steady_clock::time_point endIteration = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(endIteration - beginIteration).count() << "[s]" << std::endl;

    // std::cout << "writing to files..." << std::endl;

    // std::ofstream rhoCSV, uCSV, pCSV, tCSV; // create 4 file writers for the output
    // rhoCSV.open("rho.csv");
    // uCSV.open("u.csv");
    // pCSV.open("p.csv");
    // tCSV.open("t.csv");
    // for (int i = 0; i < iterations; i++) // loop through the iterations and points and write to file, this is very slow and should be changed
    // {
    //     for (int point = 0; point < domain.points.size(); point++)
    //     {
    //         rhoCSV << rhos[i][point];
    //         uCSV << us[i][point];
    //         pCSV << ps[i][point];
    //         if (point + 1 < domain.points.size())
    //         {
    //             rhoCSV << ",";
    //             uCSV << ",";
    //             pCSV << ",";
    //         }
    //     }
    //     rhoCSV << "\n"; // use n not endl as much faster
    //     uCSV << "\n";
    //     pCSV << "\n";
    //     tCSV << ts[i] << "\n";
    // }
    // rhoCSV.close();
    // uCSV.close();
    // pCSV.close();
    // tCSV.close();
    // std::chrono::steady_clock::time_point endSave = std::chrono::steady_clock::now();
    // std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(endSave - endIteration).count() << "[s]" << std::endl;
}

int main()
{
    std::cout << "running" << std::endl;
    sodTests(333);

    // domainTest();
}