#pragma once

#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace gslStat {
    template<typename Scalar, typename Derived>
    class RandomFiller {
    protected:
        const gsl_rng_type *T;
        gsl_rng *r;

    public:
        RandomFiller() {
            gsl_rng_env_setup();

            T = gsl_rng_default;
            r = gsl_rng_alloc(T);
        }

        ~RandomFiller() {
            gsl_rng_free(r);
        }

    private:
        Derived &derived() { return *static_cast<Derived *>(this); }

        const Derived &derived() const { return *static_cast<const Derived *>(this); }

    protected:
        Scalar generate() {
            return derived().generate();
        }

    public:
        template<typename SRC>
        void fill(SRC &src) {
            for (size_t i = 0; i < src.rows(); i++) {
                for (size_t j = 0; j < src.cols(); j++) {
                    src(i, j) = generate();
                }
            }
        }

        template<typename SRC, typename ... OTHER>
        void fill(SRC &src, OTHER &... other) {
            fill(src);
            fill(other...);
        };
    };

    class RNormal : public RandomFiller<double, RNormal> {
        using super_t = RandomFiller<double, RNormal>;
        double mean;
        double std;
    public:
        RNormal(const double _mean = 0, const double _std = 1) : super_t(), mean(_mean), std(_std) {}

        double generate() {
            return mean + gsl_ran_gaussian_ratio_method(super_t::r, std);
        }

    };

    /** Classes to generate one (continuous) uniformly distributed random variable */
    class RUnifCont : public RandomFiller<double, RUnifCont> {
        using super_t = RandomFiller<double, RUnifCont>;
        double lower;
        double upper;
        double width;
    public:
        RUnifCont(const double _lower = 0, const double _upper = 1) : super_t() {
            assert(_upper > _lower);
            lower = _lower;
            upper = _upper;
            width = upper - lower;
        }

        double generate() {
            return width * gsl_rng_uniform(r) - lower;
        }
    };

    /** Classes to generate one gamma distributed random variable */
    class RGamma : public RandomFiller<double, RGamma> {
        using super_t = RandomFiller<double, RGamma>;
        double a;
        double b;

    public:
        RGamma(double _a, double _b) : super_t(), a(_a), b(_b) {}

        double generate() {
            return gsl_ran_gamma(r, a, b);
        }
    };

    class RExp : public RandomFiller<double, RExp> {
        using super_t = RandomFiller<double, RExp>;
        double mu;

    public:
        RExp(const double _mu) : super_t(), mu(_mu) {}

        double generate() {
            return gsl_ran_exponential(r, mu);
        }
    };

    class RLaplace : public RandomFiller<double, RLaplace> {
        using super_t = RandomFiller<double, RLaplace>;
        double lambda;

    public:
        RLaplace(const double _lambda) : super_t(), lambda(_lambda) {}

        double generate() {
            return gsl_ran_laplace(r, lambda);
        }
    };

    class RExpPower : public RandomFiller<double, RExpPower> {
        using super_t = RandomFiller<double, RExpPower>;
        double a;
        double b;

    public:
        RExpPower(const double _a, const double _b) : super_t(), a(_a), b(_b) {}

        double generate() {
            return gsl_ran_exppow(r, a, b);
        }
    };

    class RCauchy : public RandomFiller<double, RCauchy> {
        using super_t = RandomFiller<double, RCauchy>;
        double a;

    public:
        RCauchy(const double _a) : super_t(), a(_a) {}

        double generate() {
            return gsl_ran_cauchy(r, a);
        }
    };

    class RRayLeigh : public RandomFiller<double, RRayLeigh> {
        using super_t = RandomFiller<double, RRayLeigh>;
        double sigma;

    public:
        RRayLeigh(const double _sigma) : super_t(), sigma(_sigma) {}

        double generate() {
            return gsl_ran_rayleigh(r, sigma);
        }
    };

    class RLogNormal : public RandomFiller<double, RLogNormal> {
        using super_t = RandomFiller<double, RLogNormal>;
        double zeta;
        double sigma;

    public:
        RLogNormal(const double _zeta, const double _sigma) : super_t(), zeta(_zeta), sigma(_sigma) {}

        double generate() {
            return gsl_ran_lognormal(r, zeta, sigma);
        }
    };

    class RChiSq : public RandomFiller<double, RChiSq> {
        using super_t = RandomFiller<double, RChiSq>;
        double nu;

    public:
        RChiSq(const double _nu) : super_t(), nu(_nu) {}

        double generate() {
            return gsl_ran_chisq(r, nu);
        }
    };

    class RFdist : public RandomFiller<double, RFdist> {
        using super_t = RandomFiller<double, RFdist>;
        double d1;
        double d2;

    public:
        RFdist(const double _d1, const double _d2) : super_t(), d1(_d1), d2(_d2) {}

        double generate() {
            return gsl_ran_fdist(r, d1, d2);
        }
    };

    class RTdist : public RandomFiller<double, RTdist> {
        using super_t = RandomFiller<double, RTdist>;
        double d;

    public:
        RTdist(const double _d) : super_t(), d(_d) {}

        double generate() {
            return gsl_ran_tdist(r, d);
        }
    };

    class RBeta : public RandomFiller<double, RBeta> {
        using super_t = RandomFiller<double, RBeta>;
        double a;
        double b;

    public:
        RBeta(const double _a, const double _b) : super_t(), a(_a), b(_b) {}

        double generate() {
            return gsl_ran_beta(r, a, b);
        }
    };

    class RLogic : public RandomFiller<double, RLogic> {
        using super_t = RandomFiller<double, RLogic>;
        double a;

    public:
        RLogic(const double _a) : super_t(), a(_a) {}

        double generate() {
            return gsl_ran_logistic(r, a);
        }
    };

}////end of namespace gslStat