#pragma once
#include <cstdlib>
#include "../gslstat.h"
using gslStat::MatrixXd;
using gslStat::VectorXd;

#include <vector>
using std::vector;

#include <iostream>
using std::cout;
using std::endl;

#include <Eigen/Dense>
#include <cmath>

#define SCALE 2.0
#define SMALL 0.1

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
namespace bOrd{
    struct Dim{
        size_t _nObs;
        size_t _dimBeta;
        size_t _dimGamma;

        Dim(const size_t nObs, const size_t dimBeta, const size_t dimGamma):_nObs(nObs), _dimBeta(dimBeta), _dimGamma(dimGamma){}
    };

    class Solver: public Dim{
    public:
        gslStat::RND rnd;

        vector<VectorXd> beta;
        VectorXd gamma;
        VectorXd b;
        VectorXd M;
        VectorXd tmpBeta;
        double t;
        ////All these have squares
        double deltaEpsilonInv;
        double deltaBetaInv;
        double deltaBInv;
        double deltaGammaInv;

        vector<vector<VectorXd>> betas;
        vector<VectorXd> gammas;
        vector<VectorXd> bs;
        vector<VectorXd> Ms;

        vector<double> ts;
        vector<double> deltaEpsilonInvs;
        vector<double> deltaBetaInvs;
        vector<double> deltaBInvs;
        vector<double> deltaGammaInvs;

        ////Tmp quantities for drawing beta_i
        VectorXd meanBeta;
        VectorXd y;
        MatrixXd varBeta;
        MatrixXd xTimeXtrans;

        ////Tmp quantities for drawing gamma
        VectorXd W;
        VectorXd meanGamma;
        MatrixXd varGamma;

        ////Tmp quantities for drawing t
        double likelihoodOld;

        ////Tmp quantities for drawing b
        VectorXd tmpB;


        ////Tmp quantities for drawing delta_beta
        VectorXd diffBBeta;
        Solver(const size_t nObs, const size_t dimBeta, const size_t dimGamma):Dim(nObs, dimBeta, dimGamma),rnd(), beta(nObs) {
            gamma = VectorXd::Random(_dimGamma);
            b = VectorXd::Random(_dimBeta);
            M = VectorXd::Random(_nObs);

            t = 1.0;
            deltaEpsilonInv = 0.5;
            deltaBetaInv = 0.5;
            deltaBInv = 0.5;
            deltaGammaInv = 0.5;

            for (size_t i =0;i<nObs;i++){
                tmpBeta = VectorXd::Random(_dimBeta);
                beta[i] = tmpBeta;
            }
            ////tmp  quantities for drawing beta_i
            meanBeta = VectorXd::Constant(_dimBeta, 0.0);
            y = VectorXd::Constant(_nObs, 0.0);
            varBeta = MatrixXd::Constant(_dimBeta, _dimBeta, 0.0);
            xTimeXtrans = MatrixXd::Constant(_dimBeta, _dimBeta, 0.0);

            ////tmp quantities for drawing gamma
            W = VectorXd::Constant(_nObs, 0.0);
            meanGamma = VectorXd::Constant(_dimGamma, 0.0);
            varGamma = MatrixXd::Constant(_dimGamma, _dimGamma, 0.0);

            ////Tmp quantities for drawing t
            likelihoodOld = 0.0;

            ////Tmp quantities for drawing b
            tmpB = VectorXd::Constant(_dimBeta, 0.0);

            ////Tmp quantities for drawing delta beta
            diffBBeta = VectorXd::Constant(_dimBeta, 0.0);
        }
        template<typename Vector>
        void get_init_likelihood(const Vector& r){
            likelihoodOld = get_likelihood(r, t);
        }
        void print() const {
            cout << "Gamma= " << endl << gamma << endl;
            cout << "b= " << endl << b << endl;
            cout << "M= " << endl << M << endl;

            cout << "t= " << t << endl;
            cout << "Inverse delta_epsilon is " << deltaEpsilonInv << endl;
            cout << "Inverse delta_beta is " << deltaBetaInv << endl;
            cout << "Inverse delta_b is " << deltaBInv << endl;
            cout << "Inverse delta_gamma is " << deltaGammaInv << endl;

            for (size_t i=0;i<_nObs;i++){
                cout << i << "'th beta is " << endl << beta[i] << endl;
            }

////Generating beta's
        }

        template<typename Matrix>
        void gen_beta(const Matrix& X, const Matrix& Z){
            gen_y(Z);

            for (size_t i=0;i<_nObs;i++){
                gen_beta_i(i, X);
            }
        }
        template<typename Matrix>
        void gen_beta_i(const size_t i, const Matrix& X){
                gen_var_beta_i(i, X);
                gen_mean_beta_i(i, X);

                rnd.rmvn(tmpBeta, meanBeta, varBeta);
                beta[i] = tmpBeta;
        }

        template<typename Matrix>
        void gen_y(const Matrix& Z) {
            for (size_t i=0;i<_nObs;i++){
                y(i) = M(i) - gslStat::innerProdd(Z.row(i), gamma);
            }
        }
                template<typename Matrix>
                void gen_mean_beta_i(const size_t i, const Matrix& X){
           meanBeta = varBeta*(deltaEpsilonInv*y(i)*X.row(i).transpose()+deltaBetaInv*b);
        }

        template<typename Matrix>
                void gen_var_beta_i(const size_t i, const Matrix& X){
            xTimeXtrans = gslStat::outerProdd(X.row(i), X.row(i));
            varBeta = (deltaEpsilonInv*xTimeXtrans+ deltaBetaInv*Eigen::MatrixXd::Identity(_dimBeta, _dimBeta)).inverse();
        }


        ////Generating gamma
        template<typename Matrix>
                void gen_gamma(const Matrix& X, const Matrix& Z){
                gen_w(X);
                gen_var_gamma(Z);
                gen_mean_gamma(Z);

                rnd.rmvn(gamma, meanGamma, varGamma);
        }

        template<typename Matrix>
        void gen_w( const Matrix& X){
            for (size_t i;i<_nObs;i++) {
                W(i) = M(i) - gslStat::innerProdd(X.row(i), beta[i]);
            }
        }

        template<typename Matrix>
                void gen_var_gamma(const Matrix& Z){
            varGamma = (deltaEpsilonInv*Z.transpose()*Z+deltaGammaInv*Eigen::MatrixXd::Identity(_dimGamma, _dimGamma)).inverse();
        }

        template<typename Matrix>
                void gen_mean_gamma(const Matrix& Z){
                meanGamma = deltaEpsilonInv*varGamma*Z.transpose()*W;
        }

        ////Generating M
        template<typename Vector, typename Matrix>
        void gen_m_i(const size_t i, const Vector& r, const Matrix& X, const Matrix& Z){
            double mean = gslStat::innerProdd(X.row(i), beta[i])+gslStat::innerProdd(Z.row(i), gamma);
            double sd = std::sqrt(1/deltaEpsilonInv);
            if (r(i)==0.0) {
                M(i) = rnd.rtn(-10000, 0, mean, sd);
            }
            else if (r(i)== 1.0){
                M(i) = rnd.rtn(0, t, mean, sd);
            }
            else M(i) = rnd.rtn(t, 10000, mean, sd);
        }

        template<typename Vector, typename Matrix>
                void gen_m(const Vector& r, const Matrix& X, const Matrix& Z){
            for (size_t i =0;i<_nObs;i++){
                gen_m_i(i,r, X, Z);
            }
        };

        template<typename Vector>
                double get_likelihood(const Vector& r, const double t){
            double result = 1.0;
            double sd = std::sqrt(1.0/deltaEpsilonInv);
            for (size_t i=0;i<_nObs;i++){
                if (r(i)==1.0){
                    result*= SCALE*(gsl_cdf_gaussian_P(M(i)-t, sd) - gsl_cdf_gaussian_P(M(i), sd));
                }
                else if (r(i)==2.0){
                    result*= SCALE*(1-gsl_cdf_gaussian_P(M(i)-t, sd));
                }
            }
            return result;
        }

        template<typename Vector>
                void gen_t(const Vector& r){
            double new_t = t + rnd.rn(0, SMALL);

            if (new_t > 0){
                double likelihoodNew = get_likelihood(r, new_t);
                double acceptProb = std::min(1.0, likelihoodNew/likelihoodOld);

                if (rnd.accept(acceptProb)){
                    t = new_t;
                    likelihoodOld = likelihoodNew;
                }

            }

        }

        void gen_b(){
            tmpB = VectorXd::Constant(_dimBeta, 0.0);

            for (size_t i = 0; i< _nObs;i++){
                tmpB += beta[i];
            }

            tmpB *= (deltaBetaInv/(_nObs*deltaBetaInv+deltaBInv));

            MatrixXd varB = (2.0/(_nObs*deltaBetaInv+deltaBInv))*Eigen::MatrixXd::Identity(_dimBeta, _dimBeta);

            rnd.rmvn(b, tmpB, varB);
        }

        template<typename Matrix>
        void gen_delta_epsilon(const Matrix& X, const Matrix& Z){
            double sum = 0.0;

            for (size_t i=0;i<_nObs;i++){
                double tmp = M(i) - gslStat::innerProdd(X.row(i), beta[i])-gslStat::innerProdd(Z.row(i), gamma)
                sum += tmp*tmp;
            }

            double a = (_nObs+1)/2.0;
            double b = 2.0/(sum+1);
            deltaEpsilonInv = rnd.rinvgamma(a, b);

        }

        void gen_delta_beta(){
            double sum =0;

            for (size_t i=0;i<_nObs;i++) {
                diffBBeta = beta[i]-b;
                sum += gslStat::innerProdd(diffBBeta, diffBBeta);

            }

            double a = (_nObs*_dimBeta+1)/2.0;
            double b = 2.0/(sum+1);

            deltaBetaInv = rnd.rinvgamma(a, b);
        }

        void gen_delta_b(){
            double sum = gslStat::innerProdd(b, b);

            double a = (_dimBeta+1.0)/2.0;
            double b = 2.0/(sum+1);

            deltaBInv = rnd.rinvgamma(a,b);
        }

        void gen_delta_gamma(){
            double sum = gslStat::innerProdd(gamma, gamma);

            double a = (_dimGamma+1.0)/2.0;
            double b = 2.0/(sum+1.0);

            deltaGammaInv = rnd.rinvgamma(a,b);
        }
    };
}