// ------ language="C++" file="src/test-rk4.cc"
#include "methods.hh"
#include "types.hh"

#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

using namespace pint;

// ------ begin <<harmonic-oscillator>>[0]
template <typename real_t, typename vector_t>
ODE<real_t, vector_t> harmonic_oscillator
    ( real_t omega_0
    , real_t zeta )
{
    return [=] (real_t t, vector_t const &y) {
        return vector_t
            ( y[1] 
            , -2 * zeta * omega_0 * y[1] - omega_0*omega_0 * y[0] );
    };
}
// ------ end

template <typename real_t, typename vector_t>
std::vector<vector_t> solve_iterative
    ( IterationStep<real_t, vector_t> step
    , std::vector<vector_t> const &y_0
    , std::vector<real_t> const &t
    , unsigned n )
{
    std::vector<vector_t> y = y_0;
    for (unsigned i = 0; i < n; ++i) {
        y = step(y, t);
    }
    return y;
}

int main()
{
    using real_t = double;
    using vector_t = Eigen::Vector2d;

    unsigned n = 21;
    auto ts = linspace<real_t>(0, 15.0, n);

    auto ode = harmonic_oscillator<real_t, vector_t>(1.0, 0.5);
    auto coarse = runge_kutta_4<real_t, vector_t>(ode);
    auto fine = iterate_step<real_t, vector_t>(coarse, 0.01); 
    auto y_0 = solve(coarse, vector_t(1.0, 0.0), ts);

    auto y = solve_iterative
        ( parareal(coarse, fine)
        , y_0
        , ts
        , 5 );

    for (unsigned i = 0; i < n; ++i) {
        std::cout << ts[i] << " " << y[i][0] << " " << y[i][1] << std::endl;
    }
    return EXIT_SUCCESS;
}
// ------ end
