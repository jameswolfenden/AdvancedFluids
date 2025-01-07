#include <gtest/gtest.h>
#include "solvers/RiemannSolver.hpp"
#include <cmath>

namespace fluid
{
    namespace testing
    {

        class RiemannSolverTest : public ::testing::Test
        {
        protected:
            RiemannSolver solver;
            Flux flux;

            // Helper function to set up a test case
            void runTest(double rhoL, double uL, double pL,
                         double rhoR, double uR, double pR,
                         double expectedP)
            {
                // Calculate sound speeds
                double aL = std::sqrt(G * pL / rhoL);
                double aR = std::sqrt(G * pR / rhoR);

                // Zero transverse velocities for 1D tests
                const double v = 0.0;
                const double w = 0.0;

                std::cout << "\nTest input:"
                          << "\nLeft:  rho=" << rhoL << ", u=" << uL << ", p=" << pL << ", a=" << aL
                          << "\nRight: rho=" << rhoR << ", u=" << uR << ", p=" << pR << ", a=" << aR << std::endl;

                bool result = solver.findStar(
                    rhoL, uL, v, w, aL, pL,
                    rhoR, uR, v, w, aR, pR,
                    flux);

                std::cout << "Resulting fluxes:"
                          << "\nf1 (mass flux): " << flux.f1
                          << "\nf2 (momentum flux): " << flux.f2
                          << "\nf5 (energy flux): " << flux.f5 << std::endl;

                ASSERT_TRUE(result) << "Solver failed to converge";

                double computedP = solver.getLastP();

                std::cout << "Expected p*: " << expectedP << std::endl;

                EXPECT_NEAR(computedP, expectedP, 1e-2)
                    << "Expected p*=" << expectedP << " but got " << computedP;
            }
        };

        // Test 1: Moderate pressure ratio
        TEST_F(RiemannSolverTest, Test1)
        {
            runTest(1.0, 0.0, 1.0,   // Left state
                    0.125, 0.0, 0.1, // Right state
                    0.30313);        // Expected p*
        }

        // Test 2: Left and right moving shock waves
        TEST_F(RiemannSolverTest, Test2)
        {
            runTest(1.0, -2.0, 0.4, // Left state
                    1.0, 2.0, 0.4,  // Right state
                    0.00189);       // Expected p*
        }

        // Test 3: Strong shock wave to the right
        TEST_F(RiemannSolverTest, Test3)
        {
            runTest(1.0, 0.0, 1000.0, // Left state
                    1.0, 0.0, 0.01,   // Right state
                    460.894);         // Expected p*
        }

        // Test 4: Strong shock wave to the left
        TEST_F(RiemannSolverTest, Test4)
        {
            runTest(1.0, 0.0, 0.01,  // Left state
                    1.0, 0.0, 100.0, // Right state
                    46.095);         // Expected p*
        }

        // Test 5: Nearly vacuum state
        TEST_F(RiemannSolverTest, Test5)
        {
            runTest(5.99924, 19.5975, 460.894,  // Left state
                    5.99242, -6.19633, 46.0950, // Right state
                    1691.64);                   // Expected p*
        }

    } // namespace testing
} // namespace fluid