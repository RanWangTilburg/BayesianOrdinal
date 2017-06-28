#include <gtest/gtest.h>
#include <iostream>

using std::cout;
using std::endl;

#include <Eigen/Dense>
using namespace Eigen;

#include <vector>
using std::vector;

class TestMisc : public ::testing::Test {

};

TEST_F(TestMisc, TESTCopyScalar){
    vector<double> vec;
    double a = 10.0;
    vec.push_back(a);
    a = 20.0;
    cout << vec[0] <<endl;
    cout << a << endl;
}

TEST_F(TestMisc, TESTCopyEigen){
    vector<MatrixXd> vec;

    MatrixXd mat = MatrixXd::Constant(5, 5, 1.0);
    vec.push_back(mat);

    mat =MatrixXd::Constant(2, 2, 0.0);
    cout << vec[0] << endl;
    cout << mat << endl;

}
