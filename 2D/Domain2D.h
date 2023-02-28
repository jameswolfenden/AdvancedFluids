#ifndef DOMAIN2D_H
#define DOMAIN2D_H

#include "Cell.h"
#include <vector>
#include "SolveRiemann.h"

class Domain2D{
    public:
    int xPhysical;
    int yPhysical;
    double xBox;
    double yBox;
    int xCellCount, xFaceCount, yCellCount, yFaceCount;
    std::vector<std::vector<Cell>> cells;
    std::vector<std::vector<Cell>> xFaces;
    std::vector<std::vector<Cell>> yFaces;
    std::vector<bool> ghostFaces;
    SolveRiemann rSolver;
    Domain2D(double xPhysical, double yPhysical, int xCellCount, int yCellCount, std::vector<bool> ghostFaces, double p, double rho, double u, double v, double dye);
    void updateCells(std::vector<Domain2D*> domain, double timeStep);
    void xfindFaces();
    void yFindFaces();
    void setGhostCells(std::vector<Domain2D*> domain);
    double timeStep();
    void xUpdateCells(double minT, std::vector<Domain2D*> domain);
    void yUpdateCells(double minT, std::vector<Domain2D*> domain);
};


#endif

