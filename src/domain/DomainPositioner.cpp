#include "domain/DomainPositioner.hpp"
#include "utils/Logger.hpp"
#include <iostream>

namespace fluid
{

    DomainPositioner::DomainPositioner(Domain *root_domain) : domain0(root_domain)
    {
        // Initialise the root domain's position
        domain0->xOrigin = -0.5 * domain0->boxDims * domain0->nx;
        domain0->yOrigin = -0.5 * domain0->boxDims * domain0->ny;
        domain0->zOrigin = -0.5 * domain0->boxDims * domain0->nz;

        // Mark root domain as processed
        passed.insert(domain0);

        // Find positions for all connected domains
        findOffsets(domain0);

        // Print final positions for debugging
        for (const auto &domain : passed)
        {
            Logger::debug("Domain " + std::to_string(domain->id) + " has origin (" +
                          std::to_string(domain->xOrigin) + ", " +
                          std::to_string(domain->yOrigin) + ", " +
                          std::to_string(domain->zOrigin) + ")");
        }
    }

    void DomainPositioner::findOffsets(Domain *domain)
    {
        // Loop through all sides 6 sides of the (3D) domain
        for (int side_index = 0; side_index < 6; side_index++)
        {
            Domain *neighbour = domain->sides[side_index];

            if (neighbour && !passed.contains(neighbour))
            {
                Logger::debug("Domain " + std::to_string(neighbour->id) +
                              " being found from " + std::to_string(domain->id) +
                              " with side " + std::to_string(side_index));

                // Mark as processed
                passed.insert(neighbour);

                // Set initial position same as parent
                neighbour->xOrigin = domain->xOrigin;
                neighbour->yOrigin = domain->yOrigin;
                neighbour->zOrigin = domain->zOrigin;

                // Adjust position based on which side it connects to
                switch (side_index)
                {
                case 0: // +x side
                    neighbour->xOrigin += domain->boxDims * (domain->nx - 1);
                    neighbour->xOrigin -= neighbour->boxDims;
                    break;
                case 1: // -x side
                    neighbour->xOrigin -= neighbour->boxDims * (neighbour->nx - 1);
                    neighbour->xOrigin += domain->boxDims;
                    break;
                case 2: // +y side
                    neighbour->yOrigin += domain->boxDims * (domain->ny - 1);
                    neighbour->yOrigin -= neighbour->boxDims;
                    break;
                case 3: // -y side
                    neighbour->yOrigin -= neighbour->boxDims * (neighbour->ny - 1);
                    neighbour->yOrigin += domain->boxDims;
                    break;
                case 4: // +z side
                    neighbour->zOrigin += domain->boxDims * (domain->nz - 1);
                    neighbour->zOrigin -= neighbour->boxDims;
                    break;
                case 5: // -z side
                    neighbour->zOrigin -= neighbour->boxDims * (neighbour->nz - 1);
                    neighbour->zOrigin += domain->boxDims;
                    break;
                }

                // Recursively position connected domains
                findOffsets(neighbour);
            }
        }
    }

} // namespace fluid