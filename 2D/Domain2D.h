#ifndef DOMAIN2D_H
#define DOMAIN2D_H

#include "Cell.h"
#include <vector>

class Domain2D{
    public:
    int xPhysical;
    int yPhysical;
    double xBox;
    double yBox;
    double minT;
    double elapsedTime;
    int xCellCount, xFaceCount, yCellCount, yFaceCount;
    std::vector<std::vector<Cell>> cells;
    std::vector<std::vector<Cell>> xFaces;
    std::vector<std::vector<Cell>> yFaces;
    Domain2D(double xPhysical, double yPhysical, int xCellCount, int yCellCount);
    void updateCells();
    void xfindFaces();
    void yFindFaces();
    void setGhostCells();
};


#endif

