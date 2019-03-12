// ------ language="C++" file="src/test-rk4.cc"
#include "runge-kutta.hh"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace pint;

template <typename real_t>
std::vector<real_t> linspace(real_t start, real_t end, unsigned n)
{
    std::vector<real_t> x(n);
    for (unsigned i = 0; i < n; ++i) {
        x[i] = start + (end * i - start * i) / (n - 1);
    }
    return x;
}

int main()
{
    using real_t = double;
    unsigned n = 100;
    auto ts = linspace<real_t>(1.0, 1.1, n);
    ODE<real_t, real_t> f = [] (real_t t, real_t y) {
        return tan(y) + 1;
    };

    auto y = solve(runge_kutta_4(f), 1.0, ts);

    for (unsigned i = 0; i < n; ++i) {
        std::cout << ts[i] << " " << y[i] << std::endl;
    }
    return EXIT_SUCCESS;
}
// ------ end
