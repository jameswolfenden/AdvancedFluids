#ifndef CONSERVED_STATE_HPP
#define CONSERVED_STATE_HPP

namespace fluid
{

    class ConservedState
    {
    public:
        ConservedState(double &rho, double &u, double &v, double &w, double &p, double &a);

        double aCalc(); // Calculate the speed of sound
        void fixVacuum(); // Make sure the state doesn't have negative density or pressure

        double u1() const;
        double u2() const;
        double u3() const;
        double u4() const;
        double u5() const;

        bool updateFromConservatives(double u1, double u2, double u3, double u4, double u5);

    private:
        double &rho;
        double &u;
        double &v;
        double &w;
        double &p;
        double &a;
    };

} // namespace fluid

#endif // CONSERVED_STATE_HPP