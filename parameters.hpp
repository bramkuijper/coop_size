#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

struct Parameters
{
    // strength of the offspring size survival function
    double survival_strength[2] = {0.5,2.0};

    // loci of the offspring size trait in each of the two environments
    double init_m = 1.0;
    double init_mb = 0.0;

    double init_bh = 0.0;

    double adult_survival = 0.1;

    // adult mortality loci
    double adult_mortality[2] = {0.5,0.5};

    double sigma[2] = {0.5,0.5};

    // kin competition strength
    double kin_comp = 0.3;

    double kin_comp_power = 2.0;

    // mutation rate of the m loci
    double mu_m = 0.01;
    double mu_bh = 0.01; // mutation rate for the bet-hedging locus
    double mu_mb = 0.01; // mutation rate for the plasticity
    double sdmu = 0.02; // standard deviation in mutation rates
    unsigned int npp = 5;
    unsigned int npatches = 5;

    // dispersal probability
    double d = 0.3;

    // parental resources
    double M[2] = {20,50};

    // whether clutch size is a fixed number
    bool fix_clutch_size = false;

    int fixed_clutch_size = 2;

    // minimum size that warrants survival
    double mmin = 1.0;

    // maximum clutch size 
    double nmax = 50;

    unsigned int max_time = 50000;

    // time span over which we calculate geometric mean fitness
    unsigned int t_geometric_fitness = 100;

    double baseline_survival = 0.1;

    unsigned int output_interval = 1;
    
    // the prefix of the file name
    std::string base_name = "sim_clutch_size";

    // whether all environments change at the same time
    bool spatially_homogenous = false;

    // if no juvs born on current patch, try to sample from
    // other patches. If one has to sample more than nsample_extinct
    // times then quit the simulation
    int nsample_extinct = 100;
}; // end of the parameter struct definition

#endif
