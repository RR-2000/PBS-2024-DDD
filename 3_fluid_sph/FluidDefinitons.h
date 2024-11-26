#ifndef PBS_FLUID_DEFINITIONS_H
#define PBS_FLUID_DEFINITIONS_H

#include <vector>
#include "Fluid.h"

namespace fluids {

    extern Fluid *boundary;

    extern Fluid *water;

    extern Fluid *honey;

    extern std::vector<Fluid*> all;

    extern std::vector<Fluid*> regularFluids;

}

#endif //PBS_FLUID_DEFINITIONS_H