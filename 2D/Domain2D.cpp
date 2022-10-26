#include <cmath>
#include <iostream>
#include "Domain2D.h"
#include "..\fluidConsts.h"

Domain2D::Domain2D(double xPhysical, double yPhysical, int xCellCount, int yCellCount)
{
    this->xPhysical = xPhysical;   // set physical size
    this->yPhysical = yPhysical;   // set physical size
    this->xCellCount = xCellCount; // set number of x cells
    this->yCellCount = yCellCount; // set number of y cells
    xFaceCount = xCellCount - 1;   // 1 less face than cells
    yFaceCount = yCellCount - 1;
    xBox = xPhysical / xCellCount; // work out size of each box
    yBox = yPhysical / yCellCount;

    elapsedTime = 0;
    minT = 0.0001; // small first time step to check waves speed is slow, should replace this with some lower bound estimate probably

    cells.resize(xCellCount, std::vector<Cell>(yCellCount));
    xFaces.resize(xFaceCount, std::vector<Cell>(yCellCount)); // the faces vectors are 'oversized', there are some faces with ghost cells on both sides but makes indexing easier
    yFaces.resize(xCellCount, std::vector<Cell>(yFaceCount));

    calcCells.resize(xCellCount, std::vector<bool>(yCellCount));

    for (int y = 0; y < yCellCount; y++)
    {
        for (int x = 0; x < xCellCount; x++)
        {
            calcCells[x][y] = true;
            if (x == 0 || y == 0 || x == xCellCount - 1 || y == yCellCount - 1)
            {
                calcCells[x][y] = false;
            }
            if (x < xCellCount / 2 && y < yCellCount / 2)
            {
                calcCells[x][y] = false;
            }
        }
    }

    for (int x = 1; x < (xCellCount - 1); x++) // exclude ghost cells
    {
        for (int y = 1; y < yCellCount - 1; y++)
        {
            cells[x][y].updatePrimatives(1, 1, 0, 0); // set inital conditions for entire domain
        }
    }

    for (int x = 1; x < (xCellCount / 3); x++) // exclude ghost cells
    {
        for (int y = 2 * yCellCount / 3; y < yCellCount - 1; y++)
        {
            cells[x][y].updatePrimatives(0.1, 0.125, 0, 0); // set some different inital conditions for some of the cells
        }
    }

    setGhostCells();
}

void Domain2D::xfindFaces()
{
    for (int y = 1; y < yFaceCount; y++)
    {
        for (int x = 0; x < xFaceCount; x++)
        {
            if (calcCells[x][y] || calcCells[x + 1][y])
            {
                // each side is a pointer to the cell on each side of the face
                Cell *xSides[2];
                xSides[0] = &cells[x][y];
                xSides[1] = &cells[x + 1][y];
                // std::cout << "x: " << x << " and " << x+1 << " , y: " << y << std::endl;
                xFaces[x][y].xFindStar(xSides); // find the star values for the half point
            }
        }
    }
}

void Domain2D::yFindFaces()
{
    for (int x = 1; x < xFaceCount; x++)
    {
        for (int y = 0; y < yFaceCount; y++)
        {
            if (calcCells[x][y] || calcCells[x][y + 1])
            {
                // each side is a pointer to the cell on each side of the face
                Cell *ySides[2];
                ySides[0] = &cells[x][y];
                ySides[1] = &cells[x][y + 1];
                // std::cout << "x: " << x << " , y: " << y << " and " << y + 1 << std::endl;
                yFaces[x][y].yFindStar(ySides); // find the star values for the half point
            }
        }
    }
}

void Domain2D::updateCells() // currently using XY, should change to XYYX or better for more accuracy
{
    xfindFaces();
    for (int y = 1; y < (yCellCount - 1); y++) // loop through all non-ghost cells
    {
        for (int x = 1; x < xCellCount - 1; x++)
        {
            if (calcCells[x][y])
            {
                // find conservatives from the half point fluxes
                double u1 = cells[x][y].u1() + minT / xBox * (xFaces[x - 1][y].f1() - xFaces[x][y].f1());
                double u2 = cells[x][y].u2() + minT / xBox * (xFaces[x - 1][y].f2() - xFaces[x][y].f2());
                double u3 = cells[x][y].u3() + minT / xBox * (xFaces[x - 1][y].f3() - xFaces[x][y].f3());
                double u4 = cells[x][y].u4() + minT / xBox * (xFaces[x - 1][y].f4() - xFaces[x][y].f4());
                cells[x][y].updateConservatives(u1, u2, u3, u4); // update the primatives in the points array from the conservatives found
            }
        }
    }
    setGhostCells();
    yFindFaces();
    double nextMinT = 9999;
    for (int x = 1; x < (xCellCount - 1); x++) // loop through all non-ghost cells
    {
        for (int y = 1; y < yCellCount - 1; y++)
        {
            if (calcCells[x][y])
            {
                // find conservatives from the half point fluxes
                double u1 = cells[x][y].u1() + minT / yBox * (yFaces[x][y - 1].g1() - yFaces[x][y].g1());
                double u2 = cells[x][y].u2() + minT / yBox * (yFaces[x][y - 1].g2() - yFaces[x][y].g2());
                double u3 = cells[x][y].u3() + minT / yBox * (yFaces[x][y - 1].g3() - yFaces[x][y].g3());
                double u4 = cells[x][y].u4() + minT / yBox * (yFaces[x][y - 1].g4() - yFaces[x][y].g4());
                cells[x][y].updateConservatives(u1, u2, u3, u4); // update the primatives in the points array from the conservatives found
                
                double C_clf = 0.9;
                if (C_clf * yBox / (std::abs(yFaces[x][y].u) + yFaces[x][y].a) < nextMinT) // set the deltaT for the next iteration being 0.9 * the min time step for this iteration
                    nextMinT = C_clf * yBox / (std::abs(yFaces[x][y].u) + yFaces[x][y].a);
                if (C_clf * xBox / (std::abs(xFaces[x][y].u) + xFaces[x][y].a) < nextMinT)
                    nextMinT = C_clf * xBox / (std::abs(xFaces[x][y].u) + xFaces[x][y].a);
            }
        }
    }
    setGhostCells();
    elapsedTime += minT;
    minT = nextMinT; // change minT for the next iteration
}

void Domain2D::setGhostCells() // set the ghost cells on all 4 sides, invert the velocities that collide with the wall
{
    /* for (int x = 1; x < (xCellCount - 1); x++) // ghost cells
     {
         if (x < xCellCount / 2)
         {
             cells[x][2 * yCellCount / 3 - 1].updatePrimatives(cells[x][1].p, cells[x][1].rho, cells[x][1].u, -cells[x][1].v);
         }
         else
         {
             cells[x][0].updatePrimatives(cells[x][1].p, cells[x][1].rho, cells[x][1].u, -cells[x][1].v);
         }
         cells[x][yCellCount - 1].updatePrimatives(cells[x][yCellCount - 2].p, cells[x][yCellCount - 2].rho, cells[x][yCellCount - 2].u, -cells[x][yCellCount - 2].v);
     }
     for (int y = 1; y < (yCellCount - 1); y++) // ghost cells
     {
         if (y < 2 * yCellCount / 3)
         {
             cells[xCellCount / 2 - 1][y].updatePrimatives(cells[1][y].p, cells[1][y].rho, -cells[1][y].u, cells[1][y].v);
         }
         else
         {
             cells[0][y].updatePrimatives(cells[1][y].p, cells[1][y].rho, -cells[1][y].u, cells[1][y].v);
         }
         cells[xCellCount - 1][y].updatePrimatives(cells[xCellCount - 2][y].p, cells[xCellCount - 2][y].rho, -cells[xCellCount - 2][y].u, cells[xCellCount - 2][y].v);
     } */

   for (int y = 0; y < (yCellCount); y++)
    {
        for (int x = 1; x < xCellCount; x++)
        {
            if (!calcCells[x - 1][y] && calcCells[x][y]) // left vertical boundary walls
            {
                cells[x-1][y].updatePrimatives(cells[x][y].p, cells[x][y].rho, -cells[x][y].u, cells[x][y].v);
            }
            if (!calcCells[x][y] && calcCells[x - 1][y]) // right vertical boundary walls
            {
                cells[x][y].updatePrimatives(cells[x-1][y].p, cells[x-1][y].rho, -cells[x-1][y].u, cells[x-1][y].v);
            }
        }
    }
    for (int x = 0; x < (xCellCount); x++)
    {
        for (int y = 1; y < yCellCount; y++)
        {
            if (!calcCells[x][y-1] && calcCells[x][y]) // 'top' vertical boundary walls
            {
                cells[x][y-1].updatePrimatives(cells[x][y].p, cells[x][y].rho, cells[x][y].u, -cells[x][y].v);
            }
            if (!calcCells[x][y] && calcCells[x][y-1]) // 'bottom' vertical boundary walls
            {
                cells[x][y].updatePrimatives(cells[x][y-1].p, cells[x][y-1].rho, cells[x][y-1].u, -cells[x][y-1].v);
            }
        }
    }
}