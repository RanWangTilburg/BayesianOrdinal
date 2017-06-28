#pragma once

#include <cassert>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../../tnormal/rtnorm.hpp"

namespace gslStat {
    class RNDBase {
    protected:
        const gsl_rng_type *T;
        gsl_rng *r;

        RNDBase() {
            gsl_rng_env_setup();

            T = gsl_rng_default;
            r = gsl_rng_alloc(T);
        }

        ~RNDBase() {
            gsl_rng_free(r);
        }

    };

    class RND : public RNDBase {
    public:

        double rn(const double mean = 0, const double sd = 1.0) {
            return mean + gsl_ran_gaussian_ratio_method(r, sd);
        }

        template<typename Matrix>
        void rn(Matrix &dst, double mean = 0, double sd = 1.0) {
            for (size_t row = 0; row < dst.rows(); row++) {
                for (size_t col = 0; col < dst.cols(); col++) {
                    dst(row, col) = rn(mean, sd);
                }
            }
        }


        template<typename Vector1, typename Vector2, typename Matrix>
        void rmvn(Vector1 &dst, const Vector2 &mean, const Matrix &var) {
            gsl_ran_multivariate_gaussian(r, mean.get_gsl_vector(), var.get_gsl_matrix(), dst.get_gsl_vector());
        };

        double rinvgamma(const double a, const double b) {
            return gsl_ran_gamma(r, a, b);
        }

        double rtn(double lower, double upper, const double mean = 0, const double sd = 1) {
            double result = rtnorm(r, lower, upper, mean, sd).first;
            return result;
        }

        double runif(double lower = 0.0, double upper = 1.0) {
            double u = (upper - lower) * gsl_rng_uniform(r) + lower;
            return u;
        }

        bool accept(double p) {
            double u = runif();

            return u < p;
        }
    };


}