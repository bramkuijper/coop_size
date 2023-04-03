
#ifndef _PATCH_HPP_
#define _PATCH_HPP_

#include <vector>
#include "individual.hpp"
#include "parameters.hpp"

class Patch
{
    public:
        bool environment_is_bad = 0;
        std::vector<Individual> breeders;
        std::vector<Individual> juveniles;

        // default constructor of a patch
        Patch(Parameters const &params);

        // copy constructor
        Patch(Patch const &other);

        void operator=(Patch const &other);
};

#endif
