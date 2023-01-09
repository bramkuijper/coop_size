#include <string>
#include "simulation.hpp"

int main(int argc, char **argv)
{
    Parameters params;

    params.sigma[0] = std::stod(argv[1]);
    params.sigma[1] = std::stod(argv[2]);
    params.M[0] = std::stod(argv[3]);
    params.M[1] = std::stod(argv[4]);
    params.mu_m = std::stod(argv[5]);
    params.mu_bh = std::stod(argv[6]);
    params.mu_mb = std::stod(argv[7]);
    params.baseline_survival = std::stod(argv[8]);
    params.survival_strength[0] = std::stod(argv[9]);
    params.survival_strength[1] = std::stod(argv[10]);
    params.max_time = std::stod(argv[11]);
    params.fix_clutch_size = std::stoi(argv[12]);
    params.kin_comp = std::stod(argv[13]);
    params.init_m = std::stod(argv[14]);
    params.base_name = argv[15];

    Simulation sim(params);

    sim.run();

    return 0;
}
