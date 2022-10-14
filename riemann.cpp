#include <iostream>
#include <cmath>
#include <fstream>
#include <errno.h>
#include "Point.h"
#include "fluidConsts.h"
#include "WavesDataPoint.h"
#include "Domain.h"
#include <vector>

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

void sodTests(int caseTest)
{
    Point left;
    Point right;
    Point* _left = &left;
    Point* _right = &right;
    if (caseTest == 1)
    {
        _left->updatePrimatives(1.0, 1.0, 0.0);
        _right->updatePrimatives(0.1, 0.125, 0);
    }
    else if (caseTest==2)
    {
        _left->updatePrimatives(0.4, 1.0, -2.0);
        _right->updatePrimatives(0.4, 1.0, 2.0);

    }
    Point *sides[2];
    sides[0] = _left;
    sides[1] = _right;

    WavesDataPoint star;
    star.findStar(sides);
    star.waveData(sides);
    saveToCSV(star, sides);

}

void domainTest()
{
    std::cout << "domaintest" << std::endl;
    Domain domain(0.1);
    int iterations = 10000;
    std::vector<std::vector<double>> rhos;
    rhos.resize(iterations, std::vector<double>(pointsSize,1));
    std::vector<std::vector<double>> us;
    us.resize(iterations, std::vector<double>(pointsSize,0));
    std::vector<std::vector<double>> ps;
    ps.resize(iterations, std::vector<double>(pointsSize,0));
    std::vector<double> ts;
    ts.resize(iterations, 0);
    std::cout << "done creating arrays" << std::endl;


    for (int i = 0; i < iterations; i++)
    {
        domain.updatePoints();
        for (int point = 0; point < domain.points.size(); point++)
        {
            rhos[i][point] = domain.points[point].rho;
            us[i][point] = domain.points[point].u;
            ps[i][point] = domain.points[point].p;
        }
        ts[i] = domain.elapsedTime;
    }

    std::ofstream rhoCSV, uCSV, pCSV, tCSV;
    rhoCSV.open("rho.csv");
    uCSV.open("u.csv");
    pCSV.open("p.csv");
    tCSV.open("t.csv");
    for (int i = 0; i < iterations; i++)
    {
        for (int point = 0; point < domain.points.size(); point++)
        {
            rhoCSV << rhos[i][point];
            uCSV << us[i][point];
            pCSV << ps[i][point];
            if (point+1<domain.points.size())
            {
            rhoCSV << ",";
            uCSV << ",";
            pCSV << ",";
            }
        }
        rhoCSV << "\n";
        uCSV << "\n";
        pCSV << "\n";
        tCSV << ts[i] << "\n";
    }
    rhoCSV.close();
    uCSV.close();
    pCSV.close();
    tCSV.close();
}

int main()
{
    //sodTests(2);
        std::cout << "running" << std::endl;


    domainTest();
}