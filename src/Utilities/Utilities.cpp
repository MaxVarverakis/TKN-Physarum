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

void Utilities::addLine(const std::string& filename, const std::vector<double>& data)
{
    std::ofstream outFile(filename + ".txt");
    checkFileOpen(outFile);

    for (unsigned int i = 0; i < data.size(); ++i)
    {
        outFile << data[i] << ',';
    }
    outFile << '\n';

    outFile.close();
}

void Utilities::addLine(const std::string& filename, const Eigen::VectorXd& data)
{
    std::ofstream outFile(filename + ".txt", std::ios::app);
    checkFileOpen(outFile);

    
    for (unsigned int i = 0; i < data.size()-1; ++i)
    {
        outFile << data(i) << ',';
    }
    outFile << data(data.size()-1);
    outFile << '\n';
}

void Utilities::parallelGraphs(std::size_t thread_ID, uint32_t seed, const float width, const float height, const unsigned int resolution, const double dt, int iter)
{
    Graph graph(seed, width, height, 2 * resolution + 1);

    // bool showCC { true };
    // bool showFC { true };
    
    while (true)
    {
        graph.evolveGraph(dt);
        if (graph.conductanceConverged() && graph.fitConverged())
        // if (graph.conductanceConverged())
        {
            // if (showCC == true)
            // {
            //     std::cout << '\n' << "Thread " << thread_ID << '\n' << "Conductance converged!" << '\n';
            //     showCC = false;
            // }
            // if (graph.fitConverged())
            // {
            //     if (showFC == true)
            //     {
            //         std::cout << '\n' << "Thread " << thread_ID << '\n' << "Fitness converged!" << '\n';
            //         showFC = false;
            //     }
                
            //     const std::string path { "/Users/max/TKN_Physarum/parallel_data_1e-4_many_graphs_" + std::to_string(resolution) + '/' };
            //     std::cout << '\n' << "Exporting CSV..." << '\n';
            //     Utilities::exportCSV(path + std::to_string(thread_ID * static_cast<size_t>(iter)), graph.sampleHSpec(100, 1e-4));
            //     std::cout << '\n' << "Export complete!" << '\n';
            //     break;
            // }
            
            std::cout << '\n' << "Thread " << thread_ID << '\n' << "Converged!" << '\n';
            
            // ************************************
            // MAKE SURE TO CREATE DIRECTORY FIRST!
            // ************************************
            
            const std::string path { "/Users/max/TKN_Physarum/parallel_data_1e-4_many_graphs_" + std::to_string(resolution) + "_clamp_5" + '/' };
            Utilities::exportCSV(path + std::to_string(thread_ID + 8 * static_cast<size_t>(iter)), graph.sampleHSpec(1000, 1e-4));
            break;
        }
    }
}

