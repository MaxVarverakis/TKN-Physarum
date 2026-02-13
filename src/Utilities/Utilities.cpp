#include "Utilities.hpp"



void Utilities::exportCSV(const std::string& filename, const std::vector<double>& data)
{
    std::ofstream outFile(filename + ".txt");
    checkFileOpen(outFile);

    for (unsigned int i = 0; i < data.size(); ++i)
    {
        outFile << data[i] << '\n';
    }

    outFile.close();
}
