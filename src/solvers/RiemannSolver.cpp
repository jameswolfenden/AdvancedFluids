#include "solvers/RiemannSolver.hpp"
#include <cmath>
#include <iostream>

namespace fluid
{
    class RiemannSolver::Impl
    {
    public:
        static void sonicRarefaction(const StateView &s, const bool &left, Flux &fl)
        {
            const int sign = left ? 1 : -1;
            const double toPow = G5 + sign * G6 * s.u / s.a;
            fl.updateFromPrimatives(s.rho * pow(toPow, G4),
                                    G5 * (sign * s.a + G7 * s.u),
                                    s.v * pow(toPow, G4),
                                    s.w * pow(toPow, G4),
                                    s.p * pow(toPow, G3));
        }

        static void solveVacuumLeft(const StateView &left, const StateView &right, Flux &fl)
        {
            if (0 >= (right.u + right.a))
            {
                std::cout << "W_R" << std::endl;
                fl.updateFromPrimatives(right);
            }
            else if (0 <= (right.u - right.a * G4))
            {
                std::cout << "W_L" << std::endl;
                fl.updateFromPrimatives(left);
            }
            else
            {
                std::cout << "W_RFan" << std::endl;
                sonicRarefaction(right, false, fl);
            }
        }

        static void solveVacuumRight(const StateView &left, const StateView &right, Flux &fl)
        {
            if (0 <= (left.u - left.a))
            {
                std::cout << "W_L" << std::endl;
                fl.updateFromPrimatives(left);
            }
            else if (0 >= (left.u + left.a * G4))
            {
                std::cout << "W_R" << std::endl;
                fl.updateFromPrimatives(right);
            }
            else
            {
                std::cout << "W_LFan" << std::endl;
                sonicRarefaction(left, true, fl);
            }
        }

        static bool checkVacuumGenerated(const StateView &left, const StateView &right, Flux &fl)
        {
            double sR = right.u - right.a * G4;
            double sL = left.u + left.a * G4;
            if (sL > sR)
            {
                return false;
            }
            if (sR > 0 && sL < 0)
            {
                std::cout << "W_0" << std::endl;
                fl.updateFromPrimatives(0.0, 0.5 * (left.u + right.u), 0.5 * (left.v + right.v), 0.5 * (left.w + right.w), 0.0);
                return true;
            }
            else if (sL >= 0)
            {
                solveVacuumLeft(left, right, fl);
                return true;
            }
            else if (sR <= 0)
            {
                solveVacuumRight(left, right, fl);
                return true;
            }
            return false;
        }

        static double calcPVRS(const StateView &left, const StateView &right)
        {
            return 0.5 * (left.p + right.p) + 0.5 * (left.u - right.u) * 0.25 * (left.rho + right.rho) * (left.a + right.a);
        }

        static double calcTwoRarefaction(const StateView &left, const StateView &right)
        {
            double p_q = pow(left.p / right.p, G1);
            double u_m = (p_q * left.u / left.a + right.u / right.a + G4 * (p_q - 1.0)) / (p_q / left.a + 1 / right.a);
            double p_TL = 1 + G7 * (left.u - u_m) / left.a;
            double p_TR = 1 + G7 * (u_m - right.u) / right.a;
            return 0.5 * (left.p * pow(p_TL, G3) + right.p * pow(p_TR, G3));
        }

        static double calcTwoShock(const StateView &left, const StateView &right, double p_PV)
        {
            double ge_L = sqrt((G5 / left.rho) / (G6 * left.p + p_PV));
            double ge_R = sqrt((G5 / right.rho) / (G6 * right.p + p_PV));
            return (ge_L * left.p + ge_R * right.p - (right.u - left.u)) / (ge_L + ge_R);
        }
    };

    bool RiemannSolver::findStar(const StateView &left, const StateView &right, Flux &fl)
    {
        if (testVacuum(left, right, fl))
        {
            std::cout << "vacuum" << std::endl;
            return true; // a vacuum is generated, values found without iteration required
        }

        double tempP;
        double tempU;
        if (!iterateP(left, right, tempP, tempU))
        {
            std::cout << "iteration failed" << std::endl;
            return false; // iteration failed, abort
        }

        if (!pickSide(left, right, fl, tempP, tempU))
        {
            std::cout << "pick side failed, rho is nan" << std::endl;
            return false; // rho is nan
        }
        lastP = tempP;
        return true;
    }

    bool RiemannSolver::pickSide(const StateView &left, const StateView &right,
                                 Flux &fl, double &tempP, double &tempU)
    {
        if (tempU >= 0.0)
        { // pick left side
            if (tempP > left.p)
                fl.updateFromPrimatives(left.rho * (((tempP / left.p) + G6) / (G6 * (tempP / left.p) + 1)),
                                        tempU, left.v, left.w, tempP);
            else
                fl.updateFromPrimatives(left.rho * pow((tempP / left.p), 1 / G), tempU, left.v, left.w, tempP);
        }
        else
        { // pick right side
            if (tempP > right.p)
                fl.updateFromPrimatives(right.rho * (((tempP / right.p) + G6) / (G6 * (tempP / right.p) + 1)),
                                        tempU, right.v, right.w, tempP);
            else
                fl.updateFromPrimatives(right.rho * pow((tempP / right.p), 1 / G), tempU, right.v, right.w, tempP);
        }
        return true;
    }

    bool RiemannSolver::testVacuum(const StateView &left, const StateView &right, Flux &fl)
    {

        bool leftVacuum = left.p == 0.0;
        bool rightVacuum = right.p == 0.0;

        if (leftVacuum && rightVacuum)
        {
            // Vacuum on both sides
            fl.updateFromPrimatives(0.0, 0.5 * (left.u + right.u), 0.5 * (left.v + right.v), 0.5 * (left.w + right.w), 0.0);
            return true;
        }
        else if (leftVacuum)
        {
            // Vacuum on left side
            Impl::solveVacuumLeft(left, right, fl);
            return true;
        }
        else if (rightVacuum)
        {
            // Vacuum on right side
            Impl::solveVacuumRight(left, right, fl);
            return true;
        }
        else if (Impl::checkVacuumGenerated(left, right, fl))
        {
            return true;
        }
        return false;
    }

    bool RiemannSolver::pickStartVal(const int errorStage, const StateView &left, const StateView &right, double &tempP)
    {
        // Implementation of pickStartVal...

        switch (errorStage)
        {
        case 0:
        {
            double p_PV = std::max(TOL, Impl::calcPVRS(left, right));
            double minP = std::min(left.p, right.p);
            double maxP = std::max(left.p, right.p);
            double maxPressureRatio = maxP / minP;
            if (maxPressureRatio < 2.0 && (minP < p_PV && p_PV < maxP))
            {
                // Select PVRS Riemann solver
                tempP = p_PV;
            }
            else if (p_PV < minP)
            {
                // Select Two-Rarefaction Riemann solver
                tempP = Impl::calcTwoRarefaction(left, right);
            }
            else
            {
                // Select Two-Shock Riemann solver with
                // PVRS as estimate
                tempP = Impl::calcTwoShock(left, right, p_PV);
            }
            break;
        }
        case 1:
        {
            tempP = 0.5 * (left.p + right.p);
            break;
        }
        case 2:
        {
            tempP = std::max(TOL, Impl::calcPVRS(left, right));
            break;
        }
        case 3:
        {
            // Select Two-Rarefaction Riemann solver
            tempP = Impl::calcTwoRarefaction(left, right);
            break;
        }
        case 4:
        {
            // Select Two-Shock Riemann solver with
            // PVRS as estimate
            double p_PV = std::max(0.0, Impl::calcPVRS(left, right));
            tempP = Impl::calcTwoShock(left, right, p_PV);
            break;
        }
        case 5:
        {
            tempP = TOL;
            break;
        }
        default:
        {
            std::cout << "Error converging on p in x" << std::endl;
            std::cout << "left.p: " << left.p << ", left.rho: " << left.rho << ", left.u: " << left.u << ", left.a: " << left.a << ", right.p: " << right.p << ", right.rho: " << right.rho << ", right.u: " << right.u << ", right.a: " << right.a << std::endl;
            return false; // the timestep is probably too small
        }
        }
        return true;
    }

    bool RiemannSolver::iterateP(const StateView &left, const StateView &right,
                                 double &tempP, double &tempU)
    {

        const int MAX_ITER = 10000;

        struct WaveState
        {
            double fs;
            double d_fs;

            static WaveState calcShockWave(const double &tempP, const double &p, const double &rho)
            {
                double A = G5 / rho;
                double B = G6 * p;
                double sqrtTerm = sqrt(A / (B + tempP));
                return {(tempP - p) * sqrtTerm, sqrtTerm * (1 - (tempP - p) / (2 * (B + tempP)))};
            }

            static WaveState calcExpansionWave(const double &tempP, const double &p, const double &rho, const double &a)
            {
                return {a * G4 * (pow(tempP / p, G1) - 1), pow(tempP / p, -G2) / (rho * a)};
            }
        };

        // Test different initial p values
        for (int errorStage = 0; errorStage < 6; errorStage++)
        {
            if (!pickStartVal(errorStage, left, right, tempP))
            {
                return false;
            }

            // Newton-Raphson iteration
            for (int iter = 0; iter < MAX_ITER; iter++)
            {
                WaveState leftWave = tempP > left.p ? WaveState::calcShockWave(tempP, left.p, left.rho) : WaveState::calcExpansionWave(tempP, left.p, left.rho, left.a);
                WaveState rightWave = tempP > right.p ? WaveState::calcShockWave(tempP, right.p, right.rho) : WaveState::calcExpansionWave(tempP, right.p, right.rho, right.a);

                double f = leftWave.fs + rightWave.fs - left.u + right.u;
                double d_f = leftWave.d_fs + rightWave.d_fs;
                double change = f / d_f;

                tempP -= change;

                // Test for convergence
                if (TOL >= 2 * fabs(change / (change + 2 * tempP)))
                {
                    tempU = 0.5 * (left.u + right.u) + 0.5 * (rightWave.fs - leftWave.fs);
                    return true;
                }
            }
        }
        return false;
    }

} // namespace fluid