#include <iostream>
#include <gtest/gtest.h>
#include <cilk/cilk.h>

class TEST_CILK : public ::testing::Test {

};

__attribute__((vector(uniform(a)))) int add(int a){
    return a+1;
}

TEST_F(TEST_CILK, TEST_CILK_FOR) {
    double data[10];

    cilk_for(int
    i = 0;
    i < 10;
    i++){
        data[i] = 0.1;
    }

    for (int i = 0; i < 10; i++) {
        EXPECT_EQ(data[i], 0.1);
    }
}

TEST_F(TEST_CILK, TEST_ARRAY_NOTATION) {
    double *data = new double[10];

    data[0:10] = 0.2;

//    for (int i=0;i<10;i++){
//        std::cout << data[i] << std::endl;
//    }

    for (int i = 0; i < 10; i++) {
        EXPECT_EQ(data[i], 0.2);
    }
    delete data;
}

TEST_F(TEST_CILK, TEST_ARRAY_REDUCE) {
    double *data = new double[10];

    data[0:10] = 0.2;
    double sum = __sec_reduce_add(data[0:10]);
    ASSERT_NEAR(sum, 2, 0.1);

}

TEST_F(TEST_CILK, TEST_SIMD) {
    double a[10], b[10];
    a[:] = 0.1;
    b[:] = 0.2;
    #pragma simd
    for (int i=0;i<10;i++){
        a[i] = a[i]+b[i];
    }

}

TEST_F(TEST_CILK, TEST_SIMD_FUNC){
    int a[10];
    a[:] = 0;

    a[:] = add(a[:]);

    for (int i=0;i<10;i++){
        std::cout << a[i] << std::endl;
    }

}