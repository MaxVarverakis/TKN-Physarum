#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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
}
