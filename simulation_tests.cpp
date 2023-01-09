#include <gtest/gtest.h>
#include "simulation.hpp"

TEST(BHSimTest, PopInitialized) {

    Parameters params;

    params.N = 5000;

    Simulation sim(params);

    EXPECT_TRUE(sim.pop.size() == params.N);

}
