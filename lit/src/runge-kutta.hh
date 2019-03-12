// ------ language="C++" file="src/runge-kutta.hh"
#pragma once
#include "types.hh"

namespace pint
{
    // ------ begin <<runge-kutta-4>>[0]
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
    // ------ end
}
// ------ end
