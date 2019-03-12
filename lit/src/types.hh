// ------ language="C++" file="src/types.hh"
#pragma once
#include <functional>
#include <vector>
#include <cassert>

namespace pint
{
    // ------ begin <<ode-type>>[0]
    template <typename real_t, typename vector_t>
    using Function = std::function
        < vector_t
          ( real_t
          , vector_t const & ) >;
    // ------ end
    // ------ begin <<ode-type>>[1]
    template <typename real_t, typename vector_t>
    using ODE = std::function
        < vector_t
          ( real_t
          , vector_t const & ) >;
    // ------ end
    // ------ begin <<ode-type>>[2]
    template <typename real_t, typename vector_t, typename matrix_t>
    using Jacobian = std::function
        < std::tuple<vector_t, matrix_t>
          ( real_t
          , vector_t const & ) >;
    // ------ end
    // ------ begin <<method-type>>[0]
    template <typename real_t, typename vector_t>
    using Integral = std::function
        < vector_t
          ( vector_t const &
          , real_t
          , real_t ) >;
    
    // ------ begin <<solve-function>>[0]
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
    // ------ end
    
    template <typename real_t, typename vector_t>
    using Method = std::function
        < Integral<real_t, vector_t> 
          ( ODE<real_t, vector_t> ) >;
    // ------ end
}
// ------ end
