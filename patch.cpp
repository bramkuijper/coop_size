#include "parameters.hpp"
#include "patch.hpp"

Patch::Patch(Parameters const &params) :
    breeders(params.npp,  // initialize stack of breeders on this patch
            Individual(
                params.init_m, 
                params.init_m, 
                params.init_mb, 
                params.init_bh))
{}

// copy constructor
Patch::Patch(Patch const &other) :
    breeders{other.breeders}
    ,environment_is_bad{other.environment_is_bad}
{}

// assignment operator
void Patch::operator=(Patch const &other)
{
    breeders = other.breeders;
    environment_is_bad = other.environment_is_bad;
}

