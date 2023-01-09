#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include "individual.hpp"
#include "parameters.hpp"


class Simulation
{
    private:
        Parameters parms;
        std::ofstream data_file;

        std::uniform_real_distribution<double> uniform;
        std::normal_distribution<double> normal;
        long unsigned time_step;
        std::random_device rd;
        unsigned int seed;
        std::mt19937 rng_r;

        // stats to track # survivors per time step
        unsigned int nsurvive = 0;

    public:
        std::vector <Individual> pop;
        std::vector <Individual> juveniles;
        std::vector <double> wvec;
        std::vector <double> varwvec;

        // the current environmental state
        bool envt;
        
        Simulation(Parameters const &params);

        // run the thing
        void run();

        // initialize the data file
        void initialize_data_files();

        void envt_change();

        void produce_offspring();

        void replace();

        void adult_survival();

        double size_dependent_survival(double const size);

        double survive_kin_competition(double const n);

        void write_parameters();

        void write_data();

        void calculate_fitness();
};

#endif
