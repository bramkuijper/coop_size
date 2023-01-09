#include <random>
#include "individual.hpp"


Individual::Individual()
{}

Individual::Individual(
        double const init_size, 
        double const init_m, 
        double const init_mb, 
        double const init_bh) :
    w{0.0}
    ,m{0.5 *init_m, 0.5*init_m}
    ,mb{0.5*init_mb,0.5*init_mb}
    ,bh{0.5*init_bh,0.5*init_bh}
    ,size{init_size}
{}

// copy constructor
Individual::Individual(Individual const &other)
{
    for (unsigned allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        m[allele_idx] = other.m[allele_idx];
        mb[allele_idx] = other.mb[allele_idx];
        bh[allele_idx] = other.bh[allele_idx];
    }

    size = other.size;
    w = other.w;
} 

// assignment operator
void Individual::operator=(Individual const &other)
{
    for (unsigned allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        m[allele_idx] = other.m[allele_idx];
        mb[allele_idx] = other.mb[allele_idx];
        bh[allele_idx] = other.bh[allele_idx];
    }

    size = other.size;
    w = other.w;
}

void Individual::mutate(Parameters const &params
        ,std::mt19937 &rng)
{
    std::uniform_real_distribution<double> uniform{0.0,1.0};

    std::normal_distribution<double> normal{0.0,params.sdmu};

    for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        // mutate bet-hedging loci
        if (uniform(rng) < params.mu_bh)
        {
            bh[allele_idx] = bh[allele_idx] + normal(rng);

            if (bh[allele_idx] < 0.0)
            {
                bh[allele_idx] = 0.0;
            }
        }

        if (uniform(rng) < params.mu_m)
        {
            m[allele_idx] = m[allele_idx] + normal(rng);

            if (m[allele_idx] < 0.0)
            {
                m[allele_idx] = 0.0;
            }
        }
        
        if (uniform(rng) < params.mu_mb)
        {
            mb[allele_idx] = mb[allele_idx] + normal(rng);
        }
    }
} // end Individual::mutate()
