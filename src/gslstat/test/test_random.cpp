//#include <gtest/gtest.h>
//#include <iostream>
//#include "../gslstat.h"
//
//class TestRandom : public ::testing::Test {
//
//};
//
//TEST_F(TestRandom, TestNormal) {
//    gslStat::MatrixXd X= gslStat::MatrixXd::Constant(10, 10, 5);
//    gslStat::MatrixXd Y = gslStat::MatrixXd::Constant(5,5,5);
//    gslStat::MatrixXd Z = gslStat::MatrixXd::Constant(5,5,0);
//    gslStat::RNormal rNormal(0, 1);
//
//    gslStat::RBeta rbeta(0.5, 0.5);
//    rNormal.fill(X, Y);
//    rbeta.fill(Z);
//    std::cout << X << std::endl;
//    std::cout << Y << std::endl;
//    std::cout << Z << std::endl;
//
//}

#include <gtest/gtest.h>
#include <iostream>

using std::cout;
using std::endl;

#include "../gslstat.h"

class TestRandomFiller : public ::testing::Test {

};

TEST_F(TestRandomFiller, TestInit) {
    gslStat::RND rnd;
}

TEST_F(TestRandomFiller, TestNormal) {
    gslStat::RND rnd;
    gslStat::MatrixXd X = gslStat::MatrixXd::Constant(5, 5, 0);
    rnd.rn(X, 0, 1);
    cout << X << endl;
}

TEST_F(TestRandomFiller, TestMVNormal) {
    gslStat::RND rnd;
    gslStat::VectorXd dst = gslStat::VectorXd::Constant(2, 0);

    gslStat::VectorXd mean = gslStat::VectorXd::Constant(2, 1);
    gslStat::MatrixXd var = gslStat::MatrixXd::Constant(2, 2, 1);
    var(1, 0) = 0;
    var(0, 1) = 0;

    rnd.rmvn(dst, mean, var);
    cout << dst << endl;
}

TEST_F(TestRandomFiller, TestInvGamma){
    gslStat::RND rnd;
    cout << rnd.rinvgamma(0.5, 0.5) << endl;

}

TEST_F(TestRandomFiller, TestTNormal){
    gslStat::RND rnd;
    double result = rnd.rtn(-1.0, 0.0, 0.0, 1.0);
    cout << result << endl;
}
