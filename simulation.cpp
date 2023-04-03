#include <cmath>
#include <cassert>
#include <random>
#include "patch.hpp"
#include "simulation.hpp"

// definition of the Simulation class constructor
Simulation::Simulation(Parameters const &params) :
    rd{}
    ,seed{rd()} // which list of random numbers to pick?
    ,rng_r{seed} // initialize the random number generator with a random_device
    ,uniform(0.0,1.0) // initialize a uniform (0,1) distribution
    ,parms{params} // initialize parameters
    ,envt{false} // initialize the environment
    ,data_file{params.base_name.c_str()} // initialize the data
    ,pop(params.npatches, Patch(params)) // initialize the population
{}

// run the actual simulation then
void Simulation::run()
{
    // initialize the output files to which the statistics will be written to
    initialize_data_files();

    // evolve the thing
    for (time_step = 0; 
            time_step < parms.max_time;
            ++time_step)
    {
        // have adults survive
        adult_survival();

        // survivors produce offspring
        produce_offspring();

        // environmental change
        envt_change();

        // calculate fitness
        calculate_fitness();

        // have offspring fill vacant breeding spots
        replace();

        if (time_step % parms.output_interval == 0)
        {
            write_data();
        }
    }
        
    write_parameters();
} // end Simulation::run()

void Simulation::calculate_fitness()
{
    double w = 0.0;
    double varw = 0.0;

    int nbreeders = 0;

    // loop over all patches
    for (std::vector<Patch>::iterator pop_iter = pop.begin();
            pop_iter != pop.end();
            ++pop_iter)
    {
        // loop over breeders within each patch
        for (std::vector <Individual>::iterator ind_iter = 
                pop_iter->breeders.begin();
                ind_iter != pop_iter->breeders.end();
                ++ind_iter)
        {
            // calculate fitness
            w += ind_iter->w;
            varw += ind_iter->w * pop_iter->w;

            nbreeders += pop_iter->breeders.size();
        }
    }

    double meanw = w / nbreeders();

    wvec.push_back(meanw);
    varwvec.push_back(varw/pop.size() - meanw * meanw);
}

// calculate stats and output them to a file
void Simulation::write_data()
{
    // to store mean and variances of offspring size
    double mean_m = 0.0;
    double var_m = 0.0;

    double mean_mb = 0.0;
    double var_mb = 0.0;

    double mean_w = 0.0;
    double var_w = 0.0;

    // mean and variance in the bet-hedging trait
    double mean_bh = 0.0;
    double var_bh = 0.0;
    
    double mean_size[2] = {0.0,0.0};
    double var_size[2] = {0.0,0.0};

    double m, mb, bh, size;

    int n_sample = 10;

    int total_n = 0;

    // loop over all patches
    for (std::vector <Patch>::iterator pop_iter = pop.begin();
            pop_iter != pop.end();
            ++pop_iter)
    {
        total_n += pop_iter->breeders.size();

        // loop over breeders within a patch
        for (std::vector <Individual>::iterator ind_iter = 
                pop_iter->breeders.begin();
                    ind_iter != pop_iter->breeders.end();
                    ++ind_iter)
        {
            // get traits
            m = ind_iter->m[0] + ind_iter->m[1];
            mb = ind_iter->mb[0] + ind_iter->mb[1];
            bh = ind_iter->bh[0] + ind_iter->bh[1];

            mean_m += m;
            var_m += m*m;

            mean_mb += mb;
            var_mb += mb*mb;

            mean_bh += bh;
            var_bh += bh*bh;

            for (int sample_idx = 0; sample_idx < n_sample; ++sample_idx)
            {
                for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
                {
                    size = m + mb * envt_idx + bh * normal(rng_r);

                    mean_size[envt_idx] += size;
                    var_size[envt_idx] += size * size;
                }
            }
        } // end for
    }

    data_file << time_step << ";" << envt << ";";

    for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
    {
        mean_size[envt_idx] /= total_n * n_sample;

        data_file << mean_size[envt_idx] << ";";
    
        var_size[envt_idx] = var_size[envt_idx] / (total_n * n_sample) - 
            mean_size[envt_idx] * mean_size[envt_idx];
        
        data_file << var_size[envt_idx] << ";";
    }


    // calculate the geometric mean fitness
    double sumlogw = 0.0;
    double sumw = 0.0;
    double ssw = 0.0;

    // calculate geometric mean fitness
    int idx_start = (int)wvec.size() - parms.t_geometric_fitness;

    if (idx_start < 0)
    {
        idx_start = 0;
    }


    for (int i = idx_start; i < wvec.size(); ++i)
    {
        sumlogw += log(wvec[i]);
        sumw += wvec[i];
        ssw += wvec[i] * wvec[i];
    }

    double wgeom = exp(sumlogw/((int)wvec.size() - idx_start));
    double wmean = sumw/((int)wvec.size() - idx_start);
    double varw = ssw/((int)wvec.size() - idx_start) - wmean * wmean;

    mean_m/=total_n;
    mean_mb/=total_n;
    mean_bh/=total_n;


    var_m= var_m/total_n - mean_m * mean_m;
    var_mb= var_mb/total_n - mean_mb * mean_mb;
    var_bh = var_bh/total_ - mean_bh * mean_bh;

    data_file << mean_m << ";";
    data_file << var_m << ";";

    data_file << mean_mb << ";";
    data_file << var_mb << ";";

    data_file << mean_bh << ";";
    data_file << var_bh << ";";

    data_file << (double)juveniles.size() << ";";
    data_file << (double)nsurvive/parms.N << ";";

    data_file << wmean << ";";
    data_file << varw << ";";
    data_file << wgeom << ";";

    data_file << std::endl;
} //end Simulation::write_data()

// initialize any data files
void Simulation::initialize_data_files()
{
    data_file << "generation;envt;";

    for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
    {
        data_file << "mean_size" << envt_idx << ";";
        data_file << "var_size" << envt_idx << ";";
    }

    data_file << "mean_m;";
    data_file << "var_m;";
    
    data_file << "mean_mb;";
    data_file << "var_mb;"; 

    data_file << "mean_bh;";
    data_file << "var_bh;";

    data_file << "juveniles;fraction_adults_survive;meanw;varw;geomw;" << std::endl;

} // end void Simulation::initialize_data_files()

// calculate survival probability based on size
double Simulation::size_dependent_survival(double const size)
{
    return(parms.baseline_survival +  
            (1.0 - parms.baseline_survival) * 
                (1.0 - exp(-parms.survival_strength[envt] * (size - parms.mmin)))
                );
} // end double Simulation::survival_prob

void Simulation::adult_survival()
{
    // reset the counter that tracks the number of surviving individuals
    // just for the sake of statistics
    nsurvive = 0;

    // auxiliary variable to store a random uniform number
    double rand_unif;

    for (std::vector<Patch>::iterator patch_iter = pop.begin();
            patch_iter != pop.end();
            ++patch_iter)
    {
        for (std::vector<Individual>::iterator breeder_iter = patch_iter->breeders.begin();
                breeder_iter != patch_iter->breeders.end();
                ++breeder_iter)
        {
            // draw a random number
            rand_unif = uniform(rng_r);

            // individual dies
            if (rand_unif > parms.adult_survival)
            {
                // note postfix decrement of patch_iter
                patch_iter->breeders.erase(breeder_iter--);
            }
        }

        nsurvive += patch_iter->breeders.size();
    } // end for std::vector Individual

    // population extinct??
    if (nsurvive < 1)
    {
        write_parameters();
        exit(1);
    }
} // end void Simulation::survive()

void Simulation::envt_change()
{
    bool change = false;

    if (parms.spatially_homogenous)
    {
        // look at environment of first patch
        if (uniform(rng_r) < parms.sigma[pop[0].environment_is_bad])
        {
            // change envts in all patches
            for (std::vector<Patch>::iterator patch_iter = pop.begin();
                    patch_iter != pop.end();
                    ++patch_iter)
            {
                patch_iter->environment_is_bad = !patch_iter->environment_is_bad;
            }
        }
    }
    else
    {
        // change envts in all patches
        for (std::vector<Patch>::iterator patch_iter = pop.begin();
                patch_iter != pop.end();
                ++patch_iter)
        {
            if (uniform(rng_r) < parms.sigma[patch_iter->environment_is_bad])
            {
                patch_iter->environment_is_bad = !patch_iter->environment_is_bad;
            }
        }
    }
}

void Simulation::replace()
{
    int n_gap = (int)parms.N - (int)pop.size();

    // random sample a whole bunch 
    std::uniform_int_distribution<int> juv_sampler(0, (int)juveniles.size() - 1);

    for (int gap_idx = 0; gap_idx < n_gap; ++gap_idx)
    {
        pop.push_back(juveniles[juv_sampler(rng_r)]);
    }

    assert(pop.size() == parms.N);
} // end void Simulation::replace

// survive kin competition
double Simulation::survive_kin_competition(double const n)
{
    return(1.0 - parms.kin_comp * n/parms.nmax);
}

// replace dead individuals with new ones
// it seems so easy: just make a fitness distribution
// but due to within individual variation in offspring size this is 
// not necessarily straightforward
void Simulation::produce_offspring()
{
    juveniles.clear();
    std::vector<double> current_clutch;

    double resources, resources_per_offspring;

    double m, plast, bh;

    size_t clutch;
    
    // make a fitness distribution
    for (std::vector<Individual>::iterator pop_iter = pop.begin();
            pop_iter != pop.end();
            ++pop_iter)
    {
        current_clutch.clear();
        resources = parms.M[envt];

        m = pop_iter->m[0] + pop_iter->m[1];
        plast = (pop_iter->mb[0] + pop_iter->mb[1]) * envt;
        bh = pop_iter->bh[0] + pop_iter->bh[1];

        
        while(true)
        {
            if (parms.fix_clutch_size && current_clutch.size() >= parms.fixed_clutch_size)
            {
                break;
            }

            resources_per_offspring = m + plast + bh * normal(rng_r);

            if (resources_per_offspring < 0)
            {
                resources_per_offspring = 0.0; 
            }


            // if Mremainder is too little to produce another offspring with
            // egg-based investment m_i (i.e., Mremainder < m_i)
            // another offspring will only be produced with a probability
            // Mremainder/m_i
            if (resources_per_offspring > resources)
            { 
                // some resources left to produce offspring
                if (uniform(rng_r) < resources/resources_per_offspring)
                {
                    current_clutch.push_back(resources_per_offspring);
                }

                // either way: if one last offspring is still being
                // produced or not, after this it is game over.
                break;
            }
            else
            {
                current_clutch.push_back(resources_per_offspring);

                // update resources
                resources -= resources_per_offspring;
            }
        } // end while true

        for (std::vector<double>::iterator clutch_iter = current_clutch.begin();
                clutch_iter != current_clutch.end();
                ++clutch_iter)
        {
            // individual survives kin comp
            if (uniform(rng_r) < size_dependent_survival(*clutch_iter) *
                    survive_kin_competition(current_clutch.size()))
            {
                Individual offspring(*pop_iter);
                offspring.mutate(parms, rng_r);
                offspring.size = *clutch_iter;
                offspring.w = 0;
                
                juveniles.push_back(offspring);

                pop_iter->w = pop_iter->w + 1.0;
            }
        }
    } // end for pop_iter

} // end void Simulation::produce_offspring()

void Simulation::write_parameters()
{
    data_file << std::endl << std::endl
        << "kin_comp;" << parms.kin_comp << std::endl
        << "kin_comp_power;" << parms.kin_comp_power << std::endl
        << "N;" << parms.N << std::endl
        << "M[0];" << parms.M[0] << std::endl
        << "M[1];" << parms.M[1] << std::endl
        << "init_m;" << parms.init_m << std::endl
        << "init_mb;" << parms.init_mb << std::endl
        << "init_bh;" << parms.init_bh << std::endl
        << "mu_m;" << parms.mu_m << std::endl
        << "mu_bh;" << parms.mu_bh << std::endl
        << "mu_mb;" << parms.mu_mb << std::endl
        << "sigma0;" << parms.sigma[0] << std::endl
        << "sigma1;" << parms.sigma[1] << std::endl
        << "fix_clutch_size;" << parms.fix_clutch_size << std::endl
        << "survival_strength0;" << parms.survival_strength[0] << std::endl
        << "survival_strength1;" << parms.survival_strength[1] << std::endl
        << "baseline_survival;" << parms.baseline_survival << std::endl
        << "seed;" << seed << std::endl;
}
