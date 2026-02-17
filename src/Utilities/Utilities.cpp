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

void Utilities::parallelGraphs(std::size_t thread_ID, uint32_t seed, const float width, const float height, const unsigned int resolution, const double dt)
{
    Graph graph(seed, width, height, 2 * resolution + 1);
    while (true)
    {
        graph.evolveGraph(dt);
        if (graph.conductanceConverged() && graph.fitConverged())
        {
            std::cout << '\n' << "Thread " << thread_ID << '\n' << "Conductance and fitness converged!" << '\n';
            const std::string path { "/Users/max/TKN_Physarum/parallel_data_15/" };
            Utilities::exportCSV(path + std::to_string(thread_ID), graph.sampleHSpec(1, 1e-2));
            break;
        }
    }
}
