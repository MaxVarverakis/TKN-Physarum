#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../Graph/Graph.hpp"

namespace Utilities
{
    void exportCSV(const std::string& filename, const std::vector<double>& data);

    template <typename FileStream>
    void checkFileOpen(const FileStream& file)
    {
        if (!file.is_open()) {
            throw std::ios_base::failure("Failed to open file!");
        }
    };

    void parallelGraphs(std::size_t thread_ID, uint32_t seed, const float width, const float height, const unsigned int resolution, const double dt);

}
