#include "domain/DomainPositioner.hpp"
#include <iostream>

namespace fluid
{

    DomainPositioner::DomainPositioner(Domain *root_domain) : domain0(root_domain)
    {
        // Initialize the root domain's position
        domain0->xOrigin = -0.5 * domain0->boxDims * domain0->nx;
        domain0->yOrigin = -0.5 * domain0->boxDims * domain0->ny;
        domain0->zOrigin = -0.5 * domain0->boxDims * domain0->nz;

        // Mark root domain as processed
        passed.insert(domain0);

        // Find positions for all connected domains
        findOffsets(domain0);

        // Print final positions for debugging
        for (auto &domain : passed)
        {
            std::cout << "Domain " << domain->id << " has origin ("
                      << domain->xOrigin << ", "
                      << domain->yOrigin << ", "
                      << domain->zOrigin << ")" << std::endl;
        }
    }

    void DomainPositioner::findOffsets(Domain *domain)
    {
        for (int side_index = 0; side_index < 6; side_index++)
        {
            Domain *neighbor = domain->sides[side_index];

            if (neighbor && !passed.contains(neighbor))
            {
                std::cout << "Domain " << neighbor->id
                          << " being found from " << domain->id
                          << " with side " << side_index << std::endl;

                // Mark as processed
                passed.insert(neighbor);

                // Set initial position same as parent
                neighbor->xOrigin = domain->xOrigin;
                neighbor->yOrigin = domain->yOrigin;
                neighbor->zOrigin = domain->zOrigin;

                // Adjust position based on which side it connects to
                switch (side_index)
                {
                case 0: // +x side
                    neighbor->xOrigin += domain->boxDims * (domain->nx - 1);
                    neighbor->xOrigin -= neighbor->boxDims;
                    break;
                case 1: // -x side
                    neighbor->xOrigin -= neighbor->boxDims * (neighbor->nx - 1);
                    neighbor->xOrigin += domain->boxDims;
                    break;
                case 2: // +y side
                    neighbor->yOrigin += domain->boxDims * (domain->ny - 1);
                    neighbor->yOrigin -= neighbor->boxDims;
                    break;
                case 3: // -y side
                    neighbor->yOrigin -= neighbor->boxDims * (neighbor->ny - 1);
                    neighbor->yOrigin += domain->boxDims;
                    break;
                case 4: // +z side
                    neighbor->zOrigin += domain->boxDims * (domain->nz - 1);
                    neighbor->zOrigin -= neighbor->boxDims;
                    break;
                case 5: // -z side
                    neighbor->zOrigin -= neighbor->boxDims * (neighbor->nz - 1);
                    neighbor->zOrigin += domain->boxDims;
                    break;
                }

                // Recursively position connected domains
                findOffsets(neighbor);
            }
        }
    }

} // namespace fluid