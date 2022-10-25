#include "Domain2D.h"
#include <iostream>
#include <fstream>

int main(){
    int xCellCount = 512;
    int yCellCount = 512;
    int iterations = 1000;
    Domain2D domain(1,1,xCellCount,yCellCount);

    std::ofstream pStartCells, pStartxFaces, pStartyFaces, pEndCells, pEndxFaces, pEndyFaces; // create 4 file writers for the output
    pStartCells.open("pStartCells.csv");
    pStartxFaces.open("pStartxFaces.csv");
    pStartyFaces.open("pStartyFaces.csv");

    for (int y = 0; y < (yCellCount); y++)
    {
        for (int x = 0; x < xCellCount; x++)
        {
            pStartCells << domain.cells[x][y].p << ",";
        }
        pStartCells << "\n";
    }
    for (int y = 0; y < (yCellCount); y++)
    {
        for (int x = 0; x < xCellCount-1; x++)
        {
            pStartxFaces << domain.xFaces[x][y].p << ",";
        }
        pStartxFaces << "\n";
    }
    for (int y = 0; y < (yCellCount-1); y++)
    {
        for (int x = 0; x < xCellCount; x++)
        {
            pStartyFaces << domain.yFaces[x][y].p << ",";
        }
        pStartyFaces << "\n";
    }

    pStartCells.close();
    pStartxFaces.close();
    pStartyFaces.close();


    for (int i = 0; i < iterations; i++)
    {
            //std::cout << domain.minT << std::endl;
            domain.updateCells();
    }
    
    pEndCells.open("pEndCells.csv");
    pEndxFaces.open("pEndxFaces.csv");
    pEndyFaces.open("pEndyFaces.csv");

    for (int y = 0; y < (yCellCount); y++)
    {
        for (int x = 0; x < xCellCount; x++)
        {
            pEndCells << domain.cells[x][y].p << ",";
        }
        pEndCells << "\n";
    }
    for (int y = 0; y < (yCellCount); y++)
    {
        for (int x = 0; x < xCellCount-1; x++)
        {
            pEndxFaces << domain.xFaces[x][y].p << ",";
        }
        pEndxFaces << "\n";
    }
    for (int y = 0; y < (yCellCount-1); y++)
    {
        for (int x = 0; x < xCellCount; x++)
        {
            pEndyFaces << domain.yFaces[x][y].p << ",";
        }
        pEndyFaces << "\n";
    }

    pEndCells.close();
    pEndxFaces.close();
    pEndyFaces.close();

    std::cout << domain.elapsedTime << std::endl;
}