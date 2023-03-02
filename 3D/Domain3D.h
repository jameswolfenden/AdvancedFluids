#ifndef DOMAIN3D_H
#define DOMAIN3D_H

#include "Box.h"
#include <vector>
#include "SolveRiemann.h"

class Domain3D{
    public:
    int xPhysical;
    int yPhysical;
    int zPhysical;
    double xBox;
    double yBox;
    double zBox;
    int xCellCount, xFaceCount, yCellCount, yFaceCount, zCellCount, zFaceCount;
    std::vector<std::vector<std::vector<Box>>> boxes;
    std::vector<std::vector<std::vector<Box>>> xFaces;
    std::vector<std::vector<std::vector<Box>>> yFaces;
    std::vector<std::vector<std::vector<Box>>> zFaces;
    std::vector<bool> ghostFaces;
    SolveRiemann rSolver;
    Domain3D(double xPhysical, double yPhysical, double zPhysical, int xCellCount, int yCellCount, int zCellCount, std::vector<bool> ghostFaces, double p, double rho, double u, double v, double w);
    void updateBoxes(std::vector<Domain3D*> domain, double timeStep);
    void xFindFaces();
    void yFindFaces();
    void zFindFaces();
    void setGhostBoxes(std::vector<Domain3D*> domain);
    double timeStep();
    void xUpdateBoxes(double minT, std::vector<Domain3D*> domain);
    void yUpdateBoxes(double minT, std::vector<Domain3D*> domain);
    void zUpdateBoxes(double minT, std::vector<Domain3D*> domain);
};


#endif

