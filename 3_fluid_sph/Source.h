#ifndef PBS_SOURCE_H
#define PBS_SOURCE_H

#include <vector>
#include "Fluid.h"

/// Base class for all sources.
/// A source is responsible for generating particles.
class Source {

public:

    /// The source will generate particles for this fluid
    Fluid* m_fluid;


    /// Initialize source for a given fluid
    explicit Source(Fluid *fluid) : m_fluid(fluid) {}


    /// Initializes the source.
    /// This method may generate particles
    virtual void init() {}


    /// Perform one time step.
    /// This method may generate particles
    /// returns true if there are more particles to be generated
    virtual bool advance(double t) {
        return false;
    };


    /// Human readable representation
    virtual std::string toString() {
        return "Source (" + m_fluid->m_name +")";
    }

};


#endif // PBS_SOURCE_H