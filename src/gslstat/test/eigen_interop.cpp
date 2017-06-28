#include <gtest/gtest.h>
#include "../gslstat.h"
#include <Eigen/Dense>
class TestEigenInterop : public ::testing::Test {

};

TEST_F(TestEigenInterop, TESTINIT) {
    gslStat::Matrix<double> mat=Eigen::MatrixXd::Constant(10, 10, 2);
    std::cout << mat << std::endl;
}

TEST_F(TestEigenInterop, TEST_PLUS){
    gslStat::Matrix<double> mat = Eigen::MatrixXd::Constant(10, 10, 2);
    gslStat::Matrix<double> mat2 = Eigen::MatrixXd::Constant(10, 10, 2);

}

TEST_F(TestEigenInterop, TEST_VECTOR){
    gslStat::VectorXd vec = gslStat::VectorXd::Constant(10, 1);
    std::cout << vec  <<std::endl;
}
