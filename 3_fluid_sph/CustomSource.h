#ifndef PBS_CUSTOM_SOURCE_H
#define PBS_CUSTOM_SOURCE_H

#include "Source.h"
#include "Fluid.h"
#include <Eigen/Dense>
#include <string>

/// Places particles based on a .xyz file
class CustomSource : public Source {
public:
    CustomSource(Fluid* fluid, std::string& file_path);

    void init() override;

    std::string m_file_path;

};


#endif // PBS_CUSTOM_SOURCE_H