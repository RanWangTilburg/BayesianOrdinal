#include <gtest/gtest.h>
#include "../bOrdinal/bOrdinal.h"

class TestBOrd : public ::testing::Test{
public:
    size_t nObs= 10;
    size_t nBeta = 2;
    size_t nGamma =3;

    bOrd::Solver solver;
    gslStat::MatrixXd X;
    gslStat::MatrixXd Z;
    gslStat::VectorXd r;
    TestBOrd():solver(nObs, nBeta, nGamma){
        X = gslStat::MatrixXd::Random(nObs, nBeta);
        Z = gslStat::MatrixXd::Random(nObs, nGamma);
        r = gslStat::VectorXd::Constant(nObs, 0);
        r(1) = 1.0;
        r(2) = 2.0;
    }
};

TEST_F(TestBOrd, DISABLED_TestInit){
    solver.print();
}

TEST_F(TestBOrd, DISABLED_TESTBeta){
    solver.print();
    cout << "After generating" << endl;

    solver.gen_beta(X, Z);
    solver.print();
}

TEST_F(TestBOrd, DISABLED_TESTGamma){
    solver.print();
    cout << "After generating" << endl;

    solver.gen_gamma(X, Z);
    solver.print();
}

TEST_F(TestBOrd, DISABLED_TESTM){
    solver.print();
    cout << "After generating" << endl;

    solver.gen_m(r, X, Z);
    solver.print();
}

TEST_F(TestBOrd, TESTLikelihood){
    cout << solver.get_likelihood(r, 1.0);
}