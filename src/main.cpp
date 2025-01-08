#include <iostream>
#include <vector>
#include "domain/Domain.hpp"
#include "domain/DomainPositioner.hpp"
#include "solvers/DomainEulerSolver.hpp"
#include "io/HDF5Writer.hpp"
#include "io/XDMFWriter.hpp"
#include "types/State.hpp"
#include "utils/Logger.hpp"

using namespace fluid;

// Function to setup the domain configuration
void setupDomains(std::vector<Domain> &domains)
{
    // Initial conditions
    State highPressure(1.32, 0.0, 0.0, 0.0, 104000);
    State ambient(1.268, 0.0, 0.0, 0.0, 101325);

    // Domain dimensions
    double cellDensity = 20;
    double pipex = 0.5;
    double pipey = 1.0;
    double pipez = 0.5;
    double boxx = 2.0;
    double boxy = 2.0;
    double boxz = 4.0;

    // Setup each domain
    domains[0].setup(0, pipex, pipey, pipez, cellDensity, highPressure);
    domains[1].setup(1, pipex, boxy, pipez, cellDensity, ambient);
    domains[2].setup(2, boxx / 2 - pipex / 2, boxy, pipez, cellDensity, ambient);
    domains[3].setup(3, pipex, boxy, boxz / 2 - pipez / 2, cellDensity, ambient);
    domains[4].setup(4, boxx / 2 - pipex / 2, boxy, pipez, cellDensity, ambient);
    domains[5].setup(5, pipex, boxy, boxz / 2 - pipez / 2, cellDensity, ambient);
    domains[6].setup(6, boxx / 2 - pipex / 2, boxy, boxz / 2 - pipez / 2, cellDensity, ambient);
    domains[7].setup(7, boxx / 2 - pipex / 2, boxy, boxz / 2 - pipez / 2, cellDensity, ambient);
    domains[8].setup(8, boxx / 2 - pipex / 2, boxy, boxz / 2 - pipez / 2, cellDensity, ambient);
    domains[9].setup(9, boxx / 2 - pipex / 2, boxy, boxz / 2 - pipez / 2, cellDensity, ambient);
}

// Function to connect domains
void connectDomains(std::vector<Domain> &domains)
{
    // 0 is +x side, 1 is -x side
    // 2 is +y side, 3 is -y side
    // 4 is +z side, 5 is -z side

    // Connect vertical pipe domains
    domains[1].sides[3] = &domains[0];
    domains[0].sides[2] = &domains[1];

    // Connect surrounding domains
    domains[1].sides[0] = &domains[2];
    domains[2].sides[1] = &domains[1];
    domains[1].sides[4] = &domains[3];
    domains[3].sides[5] = &domains[1];
    domains[1].sides[1] = &domains[4];
    domains[4].sides[0] = &domains[1];
    domains[1].sides[5] = &domains[5];
    domains[5].sides[4] = &domains[1];

    // Connect outer domains
    domains[6].sides[5] = &domains[2];
    domains[2].sides[4] = &domains[6];
    domains[6].sides[1] = &domains[3];
    domains[3].sides[0] = &domains[6];
    domains[7].sides[0] = &domains[3];
    domains[3].sides[1] = &domains[7];
    domains[7].sides[5] = &domains[4];
    domains[4].sides[4] = &domains[7];
    domains[8].sides[4] = &domains[4];
    domains[4].sides[5] = &domains[8];
    domains[8].sides[0] = &domains[5];
    domains[5].sides[1] = &domains[8];
    domains[9].sides[1] = &domains[5];
    domains[5].sides[0] = &domains[9];
    domains[9].sides[4] = &domains[2];
    domains[2].sides[5] = &domains[9];
}

int main()
{
    try
    {
        // Set log level
        Logger::setLevel(LogLevel::DEBUG);

        // Create domains
        std::vector<Domain> domains(10);

        // Setup domain geometry and initial conditions
        setupDomains(domains);

        // Connect domains
        connectDomains(domains);

        // Position domains in space
        DomainPositioner positioner(&domains[0]);

        // Setup solver
        DomainEulerSolver solver(0.7); // CFL number of 0.7

        // Setup time stepping
        std::vector<double> time{0.0};
        const double timeEnd = 0.01;
        int iteration = 0;

        // Setup output
        std::string filename("simulation");
        HDF5Writer writer(filename);

        // Write initial state
        writer.writeCoordinates(domains);
        writer.writeTimestep(domains, time.back(), iteration);

        // Main time stepping loop
        while (time.back() < timeEnd)
        {
            Logger::info("Starting iteration " + std::to_string(iteration) +
                         " at time " + std::to_string(time.back()));

            // Update solution
            if (!solver.updateDomains(domains))
            {
                Logger::error("Error in solver at iteration " + std::to_string(iteration));
                return 1;
            }

            // Update time and iteration count
            time.push_back(time.back() + 2 * solver.getMinTimeStep());
            iteration++;

            // Write output
            writer.writeTimestep(domains, time.back(), iteration);

            Logger::info("Completed iteration " + std::to_string(iteration) +
                         " at time " + std::to_string(time.back()));
        }

        // Write XDMF metadata
        XDMFWriter xdmfWriter(domains, filename, time);

        Logger::info("Simulation completed successfully!");
        return 0;
    }
    catch (const std::exception &e)
    {
        Logger::error("Error: " + std::string(e.what()));
        return 1;
    }
}