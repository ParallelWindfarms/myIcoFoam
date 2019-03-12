// ------ language="C++" file="src/parareal.hh"
#pragma once
#include "types.hh"

namespace pint
{
    // ------ begin <<time-decomposition>>[0]
    template <typename real_t>
    using TimeDecomposition = std::vector<std::tuple<real_t, real_t>>;
    // ------ end
    // ------ begin <<time-decomposition>>[1]
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
    // ------ end
}
// ------ end
