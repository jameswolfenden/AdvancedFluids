#include "Domain2D.h"
#include <iostream>
#include <fstream>

int main(){
    int xCellCount = 512;
    int yCellCount = 512;
    int iterations = 1000;
    Domain2D domain(1,1,xCellCount,yCellCount); // create 1x1 physical size domain

    std::ofstream pStartCells, pEndCells; // create file writers for the output

    pStartCells.open("pStartCells.csv");
    for (int y = 0; y < (yCellCount); y++) // write inital start cell pressures
    {
        for (int x = 0; x < xCellCount; x++)
        {
            pStartCells << domain.cells[x][y].p << ",";
        }
        pStartCells << "\n";
    }
    pStartCells.close();

    for (int i = 0; i < iterations; i++) // iterate updating the domain's cells
    {
            domain.updateCells();
    }
    
    pEndCells.open("pEndCells.csv");
    for (int y = 0; y < (yCellCount); y++) // write final cell pressures
    {
        for (int x = 0; x < xCellCount; x++)
        {
            pEndCells << domain.cells[x][y].p << ",";
        }
        pEndCells << "\n";
    }
    pEndCells.close();

    std::cout << domain.elapsedTime << std::endl; // output elapsed time (of simulation time) to console
}