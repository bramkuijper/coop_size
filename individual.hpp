#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"


class Individual
{
    public:
        double w = 0;
        double m[2] = {0.0,0.0};
        double mb[2] = {0.0,0.0};
        double bh[2] = {0.0,0.0};
        double size = 0.0;

    Individual();

    Individual(
            double const size_init
            ,double const m_init
            ,double const mb_init
            ,double const init_bh 
            );

    Individual(Individual const &other);

    void operator=(Individual const &other);

    void mutate(Parameters const &params
            ,std::mt19937 &rng);
};

#endif
