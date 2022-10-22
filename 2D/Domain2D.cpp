#include <cmath>
#include <iostream>
#include "Domain2D.h"
#include "..\fluidConsts.h"

Domain2D::Domain2D(double xPhysical, double yPhysical, int xCellCount, int yCellCount)
{
    this->xPhysical = xPhysical; // set physical size
    this->yPhysical = yPhysical; // set physical size
    this->xCellCount = xCellCount;
    this->yCellCount = yCellCount;
    xFaceCount = xCellCount - 1;
    yFaceCount = yCellCount - 1;

    xBox = xPhysical / xCellCount; // work out size of each box
    yBox = yPhysical / yCellCount;

    elapsedTime = 0;
    cells.resize(xCellCount, std::vector<Cell>(yCellCount));
    xFaces.resize(xFaceCount, std::vector<Cell>(yCellCount)); // the faces vectors are oversized, faces with ghost cells on both sides, makes indexing easier
    yFaces.resize(xCellCount, std::vector<Cell>(yFaceCount));

    for (int x = 1; x < (xCellCount - 1); x++) // exclude ghost cells
    {
        for (int y = 1; y < yCellCount - 1; y++)
        {
            cells[x][y].updatePrimatives(1, 1, 0, 0);
        }
    }

    for (int x = 1; x < (xCellCount/2); x++) // exclude ghost cells
    {
        for (int y = 1; y < yCellCount/2; y++)
        {
            cells[x][y].updatePrimatives(0.1, 0.125, 0, 0);
        }
    }

    std::cout << "i have no idea if the ghost cells are being set right (should both velocity be inverted?)" << std::endl;
    setGhostCells();
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

            //std::cout << "x face between " << x << ", " << y << " and "<< x+1 << ", " << y << std::endl;

            xFaces[x][y].xFindStar(xSides); // find the star values for the half point

            if (0.49 * xBox / (xFaces[x][y].u + xFaces[x][y].a) < minT) // change minT if the time step for this iteration should be smaller
                minT = 0.49 * xBox / (xFaces[x][y].u + xFaces[x][y].a);
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

            //std::cout << "y face between " << x << ", " << y << " and "<< x << ", " << y+1 << std::endl;


            yFaces[x][y].yFindStar(ySides); // find the star values for the half point

            if (0.49 * yBox / (yFaces[x][y].u + yFaces[x][y].a) < minT) // change minT if the time step for this iteration should be smaller
                minT = 0.49 * yBox / (yFaces[x][y].u + yFaces[x][y].a);
        }
    }
}

void Domain2D::updateCells()
{
    xfindFaces();
    minT=0.001;
    for (int y = 1; y < (yCellCount - 1); y++)
    {
        for (int x = 1; x < xCellCount - 1; x++)
        {
            // find conservatives from the half point fluxes
            double u1 = cells[x][y].u1() + minT / xBox * (xFaces[x - 1][y].f1() - xFaces[x][y].f1());
            double u2 = cells[x][y].u2() + minT / xBox * (xFaces[x - 1][y].f2() - xFaces[x][y].f2());
            double u3 = cells[x][y].u3() + minT / xBox * (xFaces[x - 1][y].f3() - xFaces[x][y].f3());
            double u4 = cells[x][y].u4() + minT / xBox * (xFaces[x - 1][y].f4() - xFaces[x][y].f4());
            cells[x][y].updateConservatives(u1, u2, u3, u4); // update the primatives in the points array from the conservatives found
            elapsedTime += minT;
            //std::cout << xFaces[x - 1][y].p <<" difference in p  on x face "<< xFaces[x][y].p << std::endl;
        }
    }
    setGhostCells();
    yFindFaces();
    minT=0.001;
    for (int x = 1; x < (xCellCount - 1); x++)
    {
        for (int y = 1; y < yCellCount - 1; y++)
        {
            // find conservatives from the half point fluxes
            double u1 = cells[x][y].u1() + minT / yBox * (yFaces[x][y - 1].g1() - yFaces[x][y].g1());
            double u2 = cells[x][y].u2() + minT / yBox * (yFaces[x][y - 1].g2() - yFaces[x][y].g2());
            double u3 = cells[x][y].u3() + minT / yBox * (yFaces[x][y - 1].g3() - yFaces[x][y].g3());
            double u4 = cells[x][y].u4() + minT / yBox * (yFaces[x][y - 1].g4() - yFaces[x][y].g4());
            cells[x][y].updateConservatives(u1, u2, u3, u4); // update the primatives in the points array from the conservatives found
            elapsedTime += minT;
              //          std::cout << yFaces[x][y-1].p <<" difference in p  on y face "<< yFaces[x][y].p << std::endl;
        }
    }
    setGhostCells();
}

void Domain2D::setGhostCells()
{
    for (int x = 1; x < (xCellCount - 1); x++) // ghost cells
    {
        cells[x][0].updatePrimatives(cells[x][1].p, cells[x][1].rho, -cells[x][1].u, -cells[x][1].v);
        cells[x][yCellCount - 1].updatePrimatives(cells[x][yCellCount - 2].p, cells[x][yCellCount - 2].rho, -cells[x][yCellCount - 2].u, -cells[x][yCellCount - 2].v);
    }
    for (int y = 1; y < (yCellCount - 1); y++) // ghost cells
    {
        cells[0][y].updatePrimatives(cells[1][y].p, cells[1][y].rho, -cells[1][y].u, -cells[1][y].v);
        cells[xCellCount - 1][y].updatePrimatives(cells[xCellCount - 2][y].p, cells[xCellCount - 2][y].rho, -cells[xCellCount - 2][y].u, -cells[xCellCount - 2][y].v);
    }
}