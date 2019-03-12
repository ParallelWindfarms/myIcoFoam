---
title: Parallel-in-time integration in OpenFOAM
author: Johan Hidding
tangle:
    prefix: pint
---

# Introduction

We'll implement Parallel-in-time methods in OpenFOAM. OpenFOAM is an extensive C++ code base with its own peculiar style of C++, lots of template libraries, different executables, and using a build system called `wmake`.

First we should have a naive implementation of a PinT method in plain C++ to make sure we understand what we're doing. We will implement Parareal on top of integrators available in the GNU Scientific Library.

# Prerequisites

I will use the following libraries:

* [GNU Scientific Library](http://gnu.org/s/gsl) generic scientific algorithms.
* [Eigen3](http://eigen.tuxfamily.org) templated vector and matrix types and linear algebra.
* [Threading Building Blocks](https://www.threadingbuildingblocks.org/) task based parallel C++ framework.

# Parareal

``` {.c++ file=src/types.hh}
#pragma once
#include <functional>
#include <vector>
#include <cassert>

namespace pint
{
    <<ode-type>>
    <<method-type>>
}
```

``` {.c++ file=src/parareal.hh}
#pragma once
#include "types.hh"

namespace pint
{
    <<time-decomposition>>
}
```

From Wikipedia:

> Parareal solves an initial value problem of the form
> 
> $$\dot{y}(t) = f(y(t), t), \quad y(t_0) = y_0 \quad \text{with} \quad t_0 \leq t \leq T.$$
>
> Here, the right hand side $f$ can correspond to the spatial discretization of a partial differential equation in a method of lines approach.

We can define the type of problem as follows

$$y_i: \mathbb{R} \to \mathbb{R},$$

``` {.c++ #ode-type}
template <typename real_t, typename vector_t>
using Function = std::function
    < vector_t
      ( real_t
      , vector_t const & ) >;
```

and the ODE in the form

$$\frac{dy_i(t)}{dt} = f_i(t, y_1(t), \dots, y_n(t))$$

``` {.c++ #ode-type}
template <typename real_t, typename vector_t>
using ODE = std::function
    < vector_t
      ( real_t
      , vector_t const & ) >;
```

and the Jacobian (when its needed by the method)

$$J_{ij} = \frac{\partial f_i(t, y(t))}{\partial y_j},$$

``` {.c++ #ode-type}
template <typename real_t, typename vector_t, typename matrix_t>
using Jacobian = std::function
    < std::tuple<vector_t, matrix_t>
      ( real_t
      , vector_t const & ) >;
```

> Parareal now requires a decomposition of the time interval $[t_0, T]$ into $P$ so-called time slices $[t_j, t_{j+1}]$ such that
>
> $$[t_0, T] = [t_0, t_1] \cup [t_1, t_2] \cup \ldots \cup [t_{P-1}, t_{P} ].$$
>
> Each time slice is assigned to one processing unit when parallelizing the algorithm, so that $P$ is equal to the number of processing units used for Parareal.

``` {.c++ #time-decomposition}
template <typename real_t>
using TimeDecomposition = std::vector<std::tuple<real_t, real_t>>;
```

We can generate a equidistant time decompostion given an interval $[t_{\rm init}, t_{\rm end}]$ and a number of slices $P$.

``` {.c++ #time-decomposition}
template <typename real_t>
TimeDecomposition<real_t> make_decomposition
    ( real_t t_init
    , real_t t_end
    , unsigned P )
{
    TimeDecomposition d;
    real_t dt = (t_end - t_init) / P;
    real_t t  = t_init;

    for (unsigned i = 0; i < P; ++i) {
        d.emplace_back(t, t + dt);
        t += dt;
    }
}
```

> Parareal is based on the iterative application of two methods for integration of ordinary differential equations. One, commonly labelled ${\mathcal {F}}$, should be of high accuracy and computational cost while the other, typically labelled ${\mathcal {G}}$, must be computationally cheap but can be much less accurate. Typically, some form of Runge-Kutta method is chosen for both coarse and fine integrator, where ${\mathcal {G}}$ might be of lower order and use a larger time step than ${\mathcal {F}}$. If the initial value problem stems from the discretization of a PDE, ${\mathcal {G}}$ can also use a coarser spatial discretization, but this can negatively impact convergence unless high order interpolation is used. The result of numerical integration with one of these methods over a time slice $[t_{j}, t_{j+1}]$ for some starting value $y_{j}$ given at $t_{j}$ is then written as
>
> $$y = \mathcal{G}(y_j, t_j, t_{j+1}).$$

The function $\mathcal{G}$ in this instance is created from the ODE, so we say $\mathcal{G}$ has the type `Integral`, while the solver (returning an integral) will have type `Method`.

``` {.c++ #method-type}
template <typename real_t, typename vector_t>
using Integral = std::function
    < vector_t
      ( vector_t const &
      , real_t
      , real_t ) >;

<<solve-function>>

template <typename real_t, typename vector_t>
using Method = std::function
    < Integral<real_t, vector_t> 
      ( ODE<real_t, vector_t> ) >;
```

In general we will want to solve the ODE for a range of times. The following convenience function takes an integral, a $y_0$ and a range of times $t_i$, and returns $y_i$.

``` {.c++ #solve-function}
template <typename real_t, typename vector_t>
std::vector<vector_t> solve
    ( Integral<real_t, vector_t> step
    , vector_t const &y_0
    , std::vector<real_t> const &t )
{
    assert(t.size() > 0);
    std::vector<vector_t> y(t.size());
    y[0] = y_0;
    for (unsigned i = 1; i < t.size(); ++i) {
        y[i] = step(y[i-1], t[i-1], t[i]);
    }
    return y;
}
```

> Serial time integration with the fine method would then correspond to a step-by-step computation of
>
> $$y_{j+1} = \mathcal{F}(y_j, t_j, t_{j+1}), \quad j=0, \ldots, P-1.$$
>
> Parareal instead uses the following iteration
>
> $$y_{j+1}^{k+1} = \mathcal{G}(y^{k+1}_j, t_j, t_{j+1}) + \mathcal{F}(y^k_j, t_j, t_{j+1}) - \mathcal{G}(y^k_j, t_j, t_{j+1}), \quad j=0, \ldots, P-1, \quad k=0, \ldots, K-1,$$
>
> where $k$ is the iteration counter. As the iteration converges and $y^{k+1}_j - y^k_j \to 0$, the terms from the coarse method cancel out and Parareal reproduces the solution that is obtained by the serial execution of the fine method only. It can be shown that Parareal converges after a maximum of $P$ iterations. For Parareal to provide speedup, however, it has to converge in a number of iterations significantly smaller than the number of time slices, that is $K \ll P$.
>
> In the Parareal iteration, the computationally expensive evaluation of $\mathcal{F}(y^k_j, t_j, t_{j+1})$ can be performed in parallel on $P$ processing units. By contrast, the dependency of $y^{k+1}_{j+1}$ on $\mathcal{G}(y^{k+1}_j, t_j, t_{j+1})$ means that the coarse correction has to be computed in serial order.

# Runge-Kutta Method

``` {.c++ file=src/runge-kutta.hh}
#pragma once
#include "types.hh"

namespace pint
{
    <<runge-kutta-4>>
}
```

``` {.c++ #runge-kutta-4}
template <typename real_t, typename vector_t>
Integral<real_t, vector_t> runge_kutta_4
    ( ODE<real_t, vector_t> f )
{
    return [=]
        ( vector_t const &y
        , real_t t_init
        , real_t t_end )
    {
        real_t   t  = t_init,
                 h  = t_end - t_init;
        vector_t k1 = h * f(t, y),
                 k2 = h * f(t + h/2, y + k1/2),
                 k3 = h * f(t + h/2, y + k2/2),
                 k4 = h * f(t + h, y + k3);
        return y + (k1 + 2*k2 + 2*k3 + k4) / 6;
    };
}
```

Let's solve the test problem

$$y' = tan(y) + 1,\ with y_0 = 1, t = [1, 1.1]$$

``` {.c++ file=src/test-rk4.cc}
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
```
