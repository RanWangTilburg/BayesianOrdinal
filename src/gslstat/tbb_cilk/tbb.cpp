#include <gtest/gtest.h>
#include <tbb/tbb.h>
#include <cilk/cilk.h>
#include <cstdio>
#include <tbb/flow_graph.h>

#pragma warning(disable: 588)
static const int N = 10000;
using namespace tbb::flow;

class TBBTest : public ::testing::Test {

};

struct square {
    int operator()(int v) { return v * v; }
};

struct cube {
    int operator()(int v) { return v * v * v; }
};

class sum {
    int &my_sum;
public:
    sum(int &s) : my_sum(s) {}

    int operator()(tuple<int, int> v) {
        my_sum += get<0>(v) + get<1>(v);
        return my_sum;
    }
};


class SumFoo {
    float *my_a;
public:
    float my_sum;

    void operator()(const tbb::blocked_range<size_t> &r) {
        float *a = my_a;
        float sum = my_sum;
        size_t end = r.end();
        for (size_t i = r.begin(); i != end; ++i)
            sum += a[i];
        my_sum = sum;
    }

    SumFoo(SumFoo &x, tbb::split) : my_a(x.my_a), my_sum(0) {}

    void join(const SumFoo &y) { my_sum += y.my_sum; }

    SumFoo(float a[]) :
            my_a(a), my_sum(0) {}
};


__attribute__((vector(uniform(input)))) int add_one(int input) {
    return input + 1;
}

void parallel_add_one(int *a) {
    tbb::parallel_for(0, N, [=](int index) { a[index] = add_one(a[index]); });
}

TEST_F(TBBTest, TBBTestForSimpleLoop) {
    int a[N];
    a[:]=0;

    parallel_add_one(a);

    for (int i = 0; i < N; i++) {
        EXPECT_EQ(a[i], 1);
    }
}

TEST_F(TBBTest, TEST_FLOW_GRAPH) {
    int result = 0;

    graph g;
    broadcast_node<int> input(g);
    function_node<int, int> squarer(g, unlimited, square());
    function_node<int, int> cuber(g, unlimited, cube());
    join_node<tuple<int, int>, queueing> join(g);
    function_node<tuple<int, int>, int>
            summer(g, serial, sum(result));

    make_edge(input, squarer);
    make_edge(input, cuber);
    make_edge(squarer, get<0>(join.input_ports()));
    make_edge(cuber, get<1>(join.input_ports()));
    make_edge(join, summer);

    for (int i = 1; i <= 10; ++i)
        input.try_put(i);
    g.wait_for_all();

    ASSERT_EQ(result, 3410);


}


