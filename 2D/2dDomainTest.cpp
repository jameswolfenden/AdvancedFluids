#include "Domain2D.h"
#include <iostream>
#include <fstream>

int main(){
    int xCellCount = 128;
    int yCellCount = 128;
    int iterations = 200;

    Domain2D domain(1,1,xCellCount,yCellCount); // create 1x1 physical size domain

    std::ofstream pStartCells, rhoStartCells, uStartCells, vStartCells, pEndCells, rhoEndCells, uEndCells, vEndCells; // create file writers for the output

    pStartCells.open("pStartCells.csv");
    rhoStartCells.open("rhoStartCells.csv");
    uStartCells.open("uStartCells.csv");
    vStartCells.open("vStartCells.csv");
    for (int y = 0; y < (yCellCount); y++) // write inital start cell pressures
    {
        for (int x = 0; x < xCellCount; x++)
        {
            pStartCells << domain.cells[x][y].p << ",";
            rhoStartCells << domain.cells[x][y].rho << ",";
            uStartCells << domain.cells[x][y].u << ",";
            vStartCells << domain.cells[x][y].v << ",";
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
            domain.updateCells();
    }
    
    pEndCells.open("pEndCells.csv");
    rhoEndCells.open("rhoEndCells.csv");
    uEndCells.open("uEndCells.csv");
    vEndCells.open("vEndCells.csv");
    for (int y = 0; y < (yCellCount); y++) // write final cell pressures
    {
        for (int x = 0; x < xCellCount; x++)
        {
            pEndCells << domain.cells[x][y].p << ",";
            rhoEndCells << domain.cells[x][y].rho << ",";
            uEndCells << domain.cells[x][y].u << ",";
            vEndCells << domain.cells[x][y].v << ",";
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

    std::cout << domain.elapsedTime << std::endl; // output elapsed time (of simulation time) to console
}