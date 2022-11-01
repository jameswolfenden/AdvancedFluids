#include <cmath>
#include <iostream>
#include "Domain2D.h"
#include "..\fluidConsts.h"

Domain2D::Domain2D(double xPhysical, double yPhysical, int xCellCount, int yCellCount, std::vector<bool> ghostFaces, double p, double rho, double u, double v)
{
    this->xPhysical = xPhysical;   // set physical size
    this->yPhysical = yPhysical;   // set physical size
    this->xCellCount = xCellCount; // set number of x cells
    this->yCellCount = yCellCount; // set number of y cells
    this->ghostFaces = ghostFaces;
    xFaceCount = xCellCount - 1; // 1 less face than cells
    yFaceCount = yCellCount - 1;
    xBox = xPhysical / xCellCount; // work out size of each box
    yBox = yPhysical / yCellCount;

    cells.resize(xCellCount, std::vector<Cell>(yCellCount));
    xFaces.resize(xFaceCount, std::vector<Cell>(yCellCount)); // the faces vectors are 'oversized', there are some faces with ghost cells on both sides but makes indexing easier
    yFaces.resize(xCellCount, std::vector<Cell>(yFaceCount));


    // this may not work if there is initial velocity as the ghost cells will not be inverted
    for (int x = 0; x < (xCellCount); x++) // include ghost cells
    {
        for (int y = 0; y < yCellCount; y++)
        {
            cells[x][y].updatePrimatives(p, rho, u, v); // set inital conditions for entire domain
        }
    }
    //setGhostCells();
}

void Domain2D::xfindFaces()
{
    for (int y = 1; y < yFaceCount; y++)
    {
        for (int x = 0; x < xFaceCount; x++)
        {
            // each side is a pointer to the cell on each side of the face
            Cell *xSides[2];
            xSides[0] = &cells[x][y];
            xSides[1] = &cells[x + 1][y];
            //std::cout << "x: " << x << " and " << x+1 << " , y: " << y << std::endl;
            xFaces[x][y].xFindStar(xSides); // find the star values for the half point
        }
    }
}

void Domain2D::yFindFaces()
{
    for (int x = 1; x < xFaceCount; x++)
    {
        for (int y = 0; y < yFaceCount; y++)
        {
            // each side is a pointer to the cell on each side of the face
            Cell *ySides[2];
            ySides[0] = &cells[x][y];
            ySides[1] = &cells[x][y + 1];
            //std::cout << "x: " << x << " , y: " << y << " and " << y + 1 << std::endl;
            yFaces[x][y].yFindStar(ySides); // find the star values for the half point
        }
    }
}

void Domain2D::updateCells(std::vector<Domain2D*> sideDomains, double minT) // currently using XY, should change to XYYX or better for more accuracy
{
    xUpdateCells(minT, sideDomains);
    yUpdateCells(minT, sideDomains);
    yUpdateCells(minT, sideDomains);
    xUpdateCells(minT, sideDomains);
}
void Domain2D::xUpdateCells(double minT, std::vector<Domain2D*> sideDomains)
{
    xfindFaces();
    for (int y = 1; y < (yCellCount - 1); y++) // loop through all non edge cells
    {
        for (int x = 1; x < xCellCount - 1; x++)
        {
            // find conservatives from the half point fluxes
            double u1 = cells[x][y].u1() + minT / xBox * (xFaces[x - 1][y].f1() - xFaces[x][y].f1());
            double u2 = cells[x][y].u2() + minT / xBox * (xFaces[x - 1][y].f2() - xFaces[x][y].f2());
            double u3 = cells[x][y].u3() + minT / xBox * (xFaces[x - 1][y].f3() - xFaces[x][y].f3());
            double u4 = cells[x][y].u4() + minT / xBox * (xFaces[x - 1][y].f4() - xFaces[x][y].f4());
            cells[x][y].updateConservatives(u1, u2, u3, u4); // update the primatives in the points array from the conservatives found
        }
    }
    setGhostCells(sideDomains);
}
void Domain2D::yUpdateCells(double minT, std::vector<Domain2D*> sideDomains)

    {
    yFindFaces();
    for (int x = 1; x < (xCellCount - 1); x++) // loop through all non edge cells
    {
        for (int y = 1; y < yCellCount - 1; y++)
        {
            // find conservatives from the half point fluxes
            double u1 = cells[x][y].u1() + minT / yBox * (yFaces[x][y - 1].g1() - yFaces[x][y].g1());
            double u2 = cells[x][y].u2() + minT / yBox * (yFaces[x][y - 1].g2() - yFaces[x][y].g2());
            double u3 = cells[x][y].u3() + minT / yBox * (yFaces[x][y - 1].g3() - yFaces[x][y].g3());
            double u4 = cells[x][y].u4() + minT / yBox * (yFaces[x][y - 1].g4() - yFaces[x][y].g4());
            cells[x][y].updateConservatives(u1, u2, u3, u4); // update the primatives in the points array from the conservatives found
        }
    }
    setGhostCells(sideDomains);
}


void Domain2D::setGhostCells(std::vector<Domain2D*> sideDomains) // set the ghost cells on all 4 sides, invert the velocities that collide with the wall
{
    for (int x = 1; x < (xCellCount - 1); x++) // ghost cells
    {
        if (ghostFaces[0])
        {
            cells[x][0].updatePrimatives(cells[x][1].p, cells[x][1].rho, cells[x][1].u, -cells[x][1].v);
        }
        else
        {
            cells[x][0].updatePrimatives(sideDomains[0]->cells[x][sideDomains[0]->yCellCount-2].p, sideDomains[0]->cells[x][sideDomains[0]->yCellCount-2].rho, sideDomains[0]->cells[x][sideDomains[0]->yCellCount-2].u, sideDomains[0]->cells[x][sideDomains[0]->yCellCount-2].v);
        }
        if (ghostFaces[1])
        {
            cells[x][yCellCount - 1].updatePrimatives(cells[x][yCellCount - 2].p, cells[x][yCellCount - 2].rho, cells[x][yCellCount - 2].u, -cells[x][yCellCount - 2].v);
        }
        else
        {
            cells[x][yCellCount - 1].updatePrimatives(sideDomains[1]->cells[x][1].p, sideDomains[1]->cells[x][1].rho, sideDomains[1]->cells[x][1].u, sideDomains[1]->cells[x][1].v);
        }
    }
    for (int y = 1; y < (yCellCount - 1); y++) // ghost cells
    {
        if (ghostFaces[2])
        {
            cells[0][y].updatePrimatives(cells[1][y].p, cells[1][y].rho, -cells[1][y].u, cells[1][y].v);
        }
        else
        {
            cells[0][y].updatePrimatives(sideDomains[2]->cells[sideDomains[2]->xCellCount-2][y].p, sideDomains[2]->cells[sideDomains[2]->xCellCount-2][y].rho, sideDomains[2]->cells[sideDomains[2]->xCellCount-2][y].u, sideDomains[2]->cells[sideDomains[2]->xCellCount-2][y].v);
        }
        if (ghostFaces[3])
        {
            cells[xCellCount - 1][y].updatePrimatives(cells[xCellCount - 2][y].p, cells[xCellCount - 2][y].rho, -cells[xCellCount - 2][y].u, cells[xCellCount - 2][y].v);
        }
        else
        {
            cells[xCellCount - 1][y].updatePrimatives(sideDomains[3]->cells[1][y].p, sideDomains[3]->cells[1][y].rho, sideDomains[3]->cells[1][y].u, sideDomains[3]->cells[1][y].v);
        }
    }
}

double Domain2D::timeStep()
{
    double step = 1000000; // high value so its larger then all starting, should change to some inital calc
    for (int x = 1; x < xCellCount - 1; x++) // loop through all non-ghost cells
    {
        for (int y = 1; y < yCellCount - 1; y++)
        {
            double C_clf = 0.5;
            // this needs optimisation
            if (C_clf * yBox / (std::abs(cells[x][y].u) + cells[x][y].a) < step) // set the deltaT for the next iteration being 0.9 * the min time step for this iteration
                step = C_clf * yBox / (std::abs(cells[x][y].u) + cells[x][y].a);
            if (C_clf * xBox / (std::abs(cells[x][y].u) + cells[x][y].a) < step) // set the deltaT for the next iteration being 0.9 * the min time step for this iteration
                step = C_clf * xBox / (std::abs(cells[x][y].u) + cells[x][y].a);
        }
    }
    return step;
}