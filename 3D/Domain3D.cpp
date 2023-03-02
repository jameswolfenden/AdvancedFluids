#include <cmath>
#include <iostream>
#include "Domain3D.h"
#include "../fluidConsts.h"
#include "SolveRiemann.h"

Domain3D::Domain3D(double xPhysical, double yPhysical, double zPhysical, int xCellCount, int yCellCount, int zCellCount, std::vector<bool> ghostFaces, double p, double rho, double u, double v, double w)
{
    this->xPhysical = xPhysical;   // set physical size
    this->yPhysical = yPhysical;   // set physical size
    this->zPhysical = zPhysical;   // set physical size
    this->xCellCount = xCellCount; // set number of x boxes
    this->yCellCount = yCellCount; // set number of y boxes
    this->zCellCount = zCellCount; // set number of z boxes
    this->ghostFaces = ghostFaces;
    xFaceCount = xCellCount - 1; // 1 less face than boxes
    yFaceCount = yCellCount - 1;
    zFaceCount = zCellCount - 1;
    xBox = xPhysical / xCellCount; // work out size of each box
    yBox = yPhysical / yCellCount;
    zBox = zPhysical / zCellCount;

    boxes.resize(xCellCount, std::vector<std::vector<Box>>(yCellCount, std::vector<Box>(zCellCount))); // create the boxes vectors
    xFaces.resize(xFaceCount, std::vector<std::vector<Box>>(yCellCount, std::vector<Box>(zCellCount))); // create the xFaces vectors
    yFaces.resize(xCellCount, std::vector<std::vector<Box>>(yFaceCount, std::vector<Box>(zCellCount))); // create the yFaces vectors
    zFaces.resize(xCellCount, std::vector<std::vector<Box>>(yCellCount, std::vector<Box>(zFaceCount))); // create the zFaces vectors    

    // this may not work if there is initial velocity as the ghost boxes will not be inverted
    for (int x = 0; x < (xCellCount); x++) // include ghost boxes
    {
        for (int y = 0; y < yCellCount; y++)
        {
            for (int z = 0; z < zCellCount; z++)
            {
                boxes[x][y][z].updatePrimatives(p, rho, u, v, w); // set all boxes to initial conditions
            }
        }
    }
    // setGhostBoxes();
}

void Domain3D::xFindFaces()
{
    for (int x = 0; x < xFaceCount; x++)
    {
        for (int y = 0; y < yCellCount; y++)
        {
            for (int z = 0; z < zCellCount; z++)
            {
                if(rSolver.findStar(boxes[x][y][z].rho, boxes[x][y][z].u, boxes[x][y][z].v, boxes[x][y][z].w, boxes[x][y][z].a, boxes[x][y][z].p, boxes[x+1][y][z].rho, boxes[x+1][y][z].u, boxes[x+1][y][z].v, boxes[x+1][y][z].w, boxes[x+1][y][z].a, boxes[x+1][y][z].p, xFaces[x][y][z].rho, xFaces[x][y][z].u, xFaces[x][y][z].v, xFaces[x][y][z].w, xFaces[x][y][z].a, xFaces[x][y][z].p))
                {
                    std::cout << "in x, x: " << x << ", y: " << y << ", z: " << z << ". Left values - rho: " << boxes[x][y][z].rho << ", u: " << boxes[x][y][z].u << ", v: " << boxes[x][y][z].v << ", w: " << boxes[x][y][z].w << ", a: " << boxes[x][y][z].a << ", p: " << boxes[x][y][z].p << ". Right - rho: " << boxes[x+1][y][z].rho << ", u: " << boxes[x+1][y][z].u << ", v: " << boxes[x+1][y][z].v << ", w: " << boxes[x+1][y][z].w << ", a: " << boxes[x+1][y][z].a << ", p: " << boxes[x+1][y][z].p << ". Result - rho: " << xFaces[x][y][z].rho << ", u: " << xFaces[x][y][z].u << ", v: " << xFaces[x][y][z].v << ", w: " << xFaces[x][y][z].w << ", a: " << xFaces[x][y][z].a << ", p: " << xFaces[x][y][z].p << std::endl;
                }
            }
        }
    }
}

void Domain3D::yFindFaces()
{
    for (int x = 0; x < xCellCount; x++)
    {
        for (int y = 0; y < yFaceCount; y++)
        {
            for (int z = 0; z < zCellCount; z++)
            {
                if(rSolver.findStar(boxes[x][y][z].rho, boxes[x][y][z].v, boxes[x][y][z].u, boxes[x][y][z].w, boxes[x][y][z].a, boxes[x][y][z].p, boxes[x][y+1][z].rho, boxes[x][y+1][z].v, boxes[x][y+1][z].u, boxes[x][y+1][z].w, boxes[x][y+1][z].a, boxes[x][y+1][z].p, yFaces[x][y][z].rho, yFaces[x][y][z].v, yFaces[x][y][z].u, yFaces[x][y][z].w, yFaces[x][y][z].a, yFaces[x][y][z].p))
                {
                    std::cout << "in y, x: " << x << ", y: " << y << ", z: " << z << ". Left values - rho: " << boxes[x][y][z].rho << ", u: " << boxes[x][y][z].u << ", v: " << boxes[x][y][z].v << ", w: " << boxes[x][y][z].w << ", a: " << boxes[x][y][z].a << ", p: " << boxes[x][y][z].p << ". Right - rho: " << boxes[x][y+1][z].rho << ", u: " << boxes[x][y+1][z].u << ", v: " << boxes[x][y+1][z].v << ", w: " << boxes[x][y+1][z].w << ", a: " << boxes[x][y+1][z].a << ", p: " << boxes[x][y+1][z].p << ". Result - rho: " << yFaces[x][y][z].rho << ", u: " << yFaces[x][y][z].u << ", v: " << yFaces[x][y][z].v << ", w: " << yFaces[x][y][z].w << ", a: " << yFaces[x][y][z].a << ", p: " << yFaces[x][y][z].p << std::endl;
                }
            }
        }
    }
}

void Domain3D::zFindFaces()
{
    for (int x = 0; x < xCellCount; x++)
    {
        for (int y = 0; y < yCellCount; y++)
        {
            for (int z = 0; z < zFaceCount; z++)
            {
                if(rSolver.findStar(boxes[x][y][z].rho, boxes[x][y][z].w, boxes[x][y][z].v, boxes[x][y][z].u, boxes[x][y][z].a, boxes[x][y][z].p, boxes[x][y][z+1].rho, boxes[x][y][z+1].w, boxes[x][y][z+1].v, boxes[x][y][z+1].u, boxes[x][y][z+1].a, boxes[x][y][z+1].p, zFaces[x][y][z].rho, zFaces[x][y][z].w, zFaces[x][y][z].v, zFaces[x][y][z].u, zFaces[x][y][z].a, zFaces[x][y][z].p))
                {
                    std::cout << "in z, x: " << x << ", y: " << y << ", z: " << z << ". Left values - rho: " << boxes[x][y][z].rho << ", u: " << boxes[x][y][z].u << ", v: " << boxes[x][y][z].v << ", w: " << boxes[x][y][z].w << ", a: " << boxes[x][y][z].a << ", p: " << boxes[x][y][z].p << ". Right - rho: " << boxes[x][y][z+1].rho << ", u: " << boxes[x][y][z+1].u << ", v: " << boxes[x][y][z+1].v << ", w: " << boxes[x][y][z+1].w << ", a: " << boxes[x][y][z+1].a << ", p: " << boxes[x][y][z+1].p << ". Result - rho: " << zFaces[x][y][z].rho << ", u: " << zFaces[x][y][z].u << ", v: " << zFaces[x][y][z].v << ", w: " << zFaces[x][y][z].w << ", a: " << zFaces[x][y][z].a << ", p: " << zFaces[x][y][z].p << std::endl;
                }
            }
        }
    }
}

void Domain3D::updateBoxes(std::vector<Domain3D *> sideDomains, double minT) // XYYX, could change to other method?
{
    setGhostBoxes(sideDomains);
    xUpdateBoxes(minT, sideDomains);
    setGhostBoxes(sideDomains);
    yUpdateBoxes(minT, sideDomains);
    setGhostBoxes(sideDomains);
    zUpdateBoxes(minT,sideDomains);
    setGhostBoxes(sideDomains);
    zUpdateBoxes(minT,sideDomains);
    setGhostBoxes(sideDomains);
    yUpdateBoxes(minT, sideDomains);
    setGhostBoxes(sideDomains);
    xUpdateBoxes(minT, sideDomains);
    setGhostBoxes(sideDomains); // dont think this ones needed since the ghosts arent used but just to make the results correct if checking

}
void Domain3D::xUpdateBoxes(double minT, std::vector<Domain3D *> sideDomains)
{
    xFindFaces();
    for (int y = 1; y < (yCellCount - 1); y++) // loop through all non edge boxes
    {
        for (int x = 1; x < xCellCount - 1; x++)
        {
            for (int z = 1; z < zCellCount - 1; z++)
            {
                // find conservatives from the half point fluxes
                double u1 = boxes[x][y][z].u1() + minT / xBox * (xFaces[x - 1][y][z].f1() - xFaces[x][y][z].f1());
                double u2 = boxes[x][y][z].u2() + minT / xBox * (xFaces[x - 1][y][z].f2() - xFaces[x][y][z].f2());
                double u3 = boxes[x][y][z].u3() + minT / xBox * (xFaces[x - 1][y][z].f3() - xFaces[x][y][z].f3());
                double u4 = boxes[x][y][z].u4() + minT / xBox * (xFaces[x - 1][y][z].f4() - xFaces[x][y][z].f4());
                double u5 = boxes[x][y][z].u5() + minT / xBox * (xFaces[x - 1][y][z].f5() - xFaces[x][y][z].f5());
                boxes[x][y][z].updateConservatives(u1, u2, u3, u4, u5);
            }
        }
    }
}
void Domain3D::yUpdateBoxes(double minT, std::vector<Domain3D *> sideDomains)
{
    yFindFaces();
    for (int x = 1; x < (xCellCount - 1); x++) // loop through all non edge boxes
    {
        for (int y = 1; y < yCellCount - 1; y++)
        {
            for (int z = 1; z < zCellCount - 1; z++)
            {
                // find conservatives from the half point fluxes
                double u1 = boxes[x][y][z].u1() + minT / yBox * (yFaces[x][y - 1][z].g1() - yFaces[x][y][z].g1());
                double u2 = boxes[x][y][z].u2() + minT / yBox * (yFaces[x][y - 1][z].g2() - yFaces[x][y][z].g2());
                double u3 = boxes[x][y][z].u3() + minT / yBox * (yFaces[x][y - 1][z].g3() - yFaces[x][y][z].g3());
                double u4 = boxes[x][y][z].u4() + minT / yBox * (yFaces[x][y - 1][z].g4() - yFaces[x][y][z].g4());
                double u5 = boxes[x][y][z].u5() + minT / yBox * (yFaces[x][y - 1][z].g5() - yFaces[x][y][z].g5());
                boxes[x][y][z].updateConservatives(u1, u2, u3, u4, u5);
            }
        }
    }
}
void Domain3D::zUpdateBoxes(double minT, std::vector<Domain3D *> sideDomains)
{
    zFindFaces();
    for (int x = 1; x < (xCellCount - 1); x++) // loop through all non edge boxes
    {
        for (int y = 1; y < yCellCount - 1; y++)
        {
            for (int z = 1; z < zCellCount - 1; z++)
            {
                // find conservatives from the half point fluxes
                double u1 = boxes[x][y][z].u1() + minT / zBox * (zFaces[x][y][z - 1].h1() - zFaces[x][y][z].h1());
                double u2 = boxes[x][y][z].u2() + minT / zBox * (zFaces[x][y][z - 1].h2() - zFaces[x][y][z].h2());
                double u3 = boxes[x][y][z].u3() + minT / zBox * (zFaces[x][y][z - 1].h3() - zFaces[x][y][z].h3());
                double u4 = boxes[x][y][z].u4() + minT / zBox * (zFaces[x][y][z - 1].h4() - zFaces[x][y][z].h4());
                double u5 = boxes[x][y][z].u5() + minT / zBox * (zFaces[x][y][z - 1].h5() - zFaces[x][y][z].h5());
                boxes[x][y][z].updateConservatives(u1, u2, u3, u4, u5);
            }
        }
    }
}

void Domain3D::setGhostBoxes(std::vector<Domain3D *> sideDomains) // set the ghost boxes on all 6 sides, invert the velocities that collide with the wall
{    
    for (int x = 1; x < (xCellCount - 1); x++) // ghost boxes
    {
        for (int z = 1; z < (zCellCount - 1); z++) 
        {
            if (ghostFaces[0])
                boxes[x][0][z].updatePrimatives(boxes[x][1][z].p, boxes[x][1][z].rho, boxes[x][1][z].u, -boxes[x][1][z].v, boxes[x][1][z].w);
            else
                boxes[x][0][z].updatePrimatives(sideDomains[0]->boxes[x][sideDomains[0]->yCellCount - 2][z].p, sideDomains[0]->boxes[x][sideDomains[0]->yCellCount - 2][z].rho, sideDomains[0]->boxes[x][sideDomains[0]->yCellCount - 2][z].u, sideDomains[0]->boxes[x][sideDomains[0]->yCellCount - 2][z].v, sideDomains[0]->boxes[x][sideDomains[0]->yCellCount - 2][z].w);
            if (ghostFaces[1])
                boxes[x][yCellCount - 1][z].updatePrimatives(boxes[x][yCellCount - 2][z].p, boxes[x][yCellCount - 2][z].rho, boxes[x][yCellCount - 2][z].u, -boxes[x][yCellCount - 2][z].v, boxes[x][yCellCount - 2][z].w);
            else
                boxes[x][yCellCount - 1][z].updatePrimatives(sideDomains[1]->boxes[x][1][z].p, sideDomains[1]->boxes[x][1][z].rho, sideDomains[1]->boxes[x][1][z].u, sideDomains[1]->boxes[x][1][z].v, sideDomains[1]->boxes[x][1][z].w);
        }
    }
    for (int y = 1; y < (yCellCount - 1); y++) // ghost boxes
    {
        for (int z = 1; z < (zCellCount - 1); z++) 
        {
            if (ghostFaces[2])
                boxes[0][y][z].updatePrimatives(boxes[1][y][z].p, boxes[1][y][z].rho, -boxes[1][y][z].u, boxes[1][y][z].v, boxes[1][y][z].w);
            else
                boxes[0][y][z].updatePrimatives(sideDomains[2]->boxes[sideDomains[2]->xCellCount - 2][y][z].p, sideDomains[2]->boxes[sideDomains[2]->xCellCount - 2][y][z].rho, sideDomains[2]->boxes[sideDomains[2]->xCellCount - 2][y][z].u, sideDomains[2]->boxes[sideDomains[2]->xCellCount - 2][y][z].v, sideDomains[2]->boxes[sideDomains[2]->xCellCount - 2][y][z].w);
            if (ghostFaces[3])
                boxes[xCellCount - 1][y][z].updatePrimatives(boxes[xCellCount - 2][y][z].p, boxes[xCellCount - 2][y][z].rho, -boxes[xCellCount - 2][y][z].u, boxes[xCellCount - 2][y][z].v, boxes[xCellCount - 2][y][z].w);
            else
                boxes[xCellCount - 1][y][z].updatePrimatives(sideDomains[3]->boxes[1][y][z].p, sideDomains[3]->boxes[1][y][z].rho, sideDomains[3]->boxes[1][y][z].u, sideDomains[3]->boxes[1][y][z].v, sideDomains[3]->boxes[1][y][z].w);
        }
    }
    for (int x = 1; x < (xCellCount - 1); x++) // ghost boxes
    {
        for (int y = 1; y < (yCellCount - 1); y++) 
        {
            if (ghostFaces[4])
                boxes[x][y][0].updatePrimatives(boxes[x][y][1].p, boxes[x][y][1].rho, boxes[x][y][1].u, boxes[x][y][1].v, -boxes[x][y][1].w);
            else
                boxes[x][y][0].updatePrimatives(sideDomains[4]->boxes[x][y][sideDomains[4]->zCellCount - 2].p, sideDomains[4]->boxes[x][y][sideDomains[4]->zCellCount - 2].rho, sideDomains[4]->boxes[x][y][sideDomains[4]->zCellCount - 2].u, sideDomains[4]->boxes[x][y][sideDomains[4]->zCellCount - 2].v, sideDomains[4]->boxes[x][y][sideDomains[4]->zCellCount - 2].w);
            if (ghostFaces[5])
                boxes[x][y][zCellCount - 1].updatePrimatives(boxes[x][y][zCellCount - 2].p, boxes[x][y][zCellCount - 2].rho, boxes[x][y][zCellCount - 2].u, boxes[x][y][zCellCount - 2].v, -boxes[x][y][zCellCount - 2].w);
            else
                boxes[x][y][zCellCount - 1].updatePrimatives(sideDomains[5]->boxes[x][y][1].p, sideDomains[5]->boxes[x][y][1].rho, sideDomains[5]->boxes[x][y][1].u, sideDomains[5]->boxes[x][y][1].v, sideDomains[5]->boxes[x][y][1].w);
        }
    }
}

double Domain3D::timeStep()
{
    double step = 1000000;                   // high value so its larger then all starting, should change to some inital calc
    for (int x = 1; x < xCellCount - 1; x++) // loop through all non-ghost boxes
    {
        for (int y = 1; y < yCellCount - 1; y++)
        {
            for (int z = 1; z < zCellCount - 1; z++)
            {
                double C_clf = 0.5;
                // this needs optimisation
                if (C_clf * yBox / (std::abs(boxes[x][y][z].v) + boxes[x][y][z].a) < step) // set the deltaT for the next iteration being 0.9 * the min time step for this iteration
                    step = C_clf * yBox / (std::abs(boxes[x][y][z].v) + boxes[x][y][z].a);
                if (C_clf * xBox / (std::abs(boxes[x][y][z].u) + boxes[x][y][z].a) < step) // set the deltaT for the next iteration being 0.9 * the min time step for this iteration
                    step = C_clf * xBox / (std::abs(boxes[x][y][z].u) + boxes[x][y][z].a);
                if (C_clf * zBox / (std::abs(boxes[x][y][z].w) + boxes[x][y][z].a) < step) // set the deltaT for the next iteration being 0.9 * the min time step for this iteration
                    step = C_clf * zBox / (std::abs(boxes[x][y][z].w) + boxes[x][y][z].a);
            }
        }
    }
    return step;
}