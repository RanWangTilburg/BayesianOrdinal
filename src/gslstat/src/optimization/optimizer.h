#pragma once

#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>

namespace gslStat{
    template<typename Derived>
    class MinimizerBase{
    private:
        /** \returns a reference to the derived object */
        Derived &derived() { return *static_cast<Derived *>(this); }

        /** \returns a const reference to the derived object */
        const Derived &derived() const { return *static_cast<const Derived *>(this); }

        const gsl_multimin_fdfminimizer_type *T;
        gsl_multimin_fdfminimizer *s;
    public:
        MinimizerBase(){
            T = gsl_multimin_fdfminimizer_vector_bfgs2;
            s = gsl_multimin_fdfminimizer_alloc (T, 2);
        }
    };

    class MinimizerWithDer:public MinimizerBase<MinimizerWithDer>{
    public:
        using super_t = MinimizerBase<MinimizerWithDer>;

    };

}////end of namespace gslStat