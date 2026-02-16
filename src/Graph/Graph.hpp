#pragma once

#include <vector>
#include <glm/glm.hpp>
#include <Sparse>
#include <random>
#include <iostream>
#include <algorithm>

struct Edge
{
    unsigned int i, j; // node indices
    
    // only stores initial values of D & Q
    // this can be made better if create Dvec in place of initializing edge.D and then populating Dvec afterwards (in regularLattice function)
    double D; // conductance
    double Q; // flow
};

struct Node
{
    glm::fvec2 pos; // for rendering circles/lines
    // s and p stored in `Eigen` vectors
};

class Graph
{
private:
    // Randomness
    uint32_t m_master_seed;
    std::mt19937 m_rng_sources;
    std::mt19937 m_rng_initD;

    unsigned int m_resolution;
    unsigned int m_sink_idx;
    std::vector<unsigned int> m_source_ids;
    const unsigned int n_sources { 30 };
    // const unsigned int n_active_node_pairs { 6 };
    const double I0;
    const double D0 { 0.001 };
    bool solverInitialized { false };
    bool fitnessConverged { false };
    const double m_tol { 1e-6 }; // for convergence
    const double c_t { 2.0 };

    std::vector<Node> m_nodes;
    std::vector<Edge> m_edges;
    
    Eigen::SparseMatrix<double> L; // Laplacian
    Eigen::VectorXd p, s; // pressures, sources/sinks
    
    // reduced versions
    Eigen::SparseMatrix<double> Lr;
    Eigen::VectorXd sr;
    
    Eigen::VectorXd Dvec, Qvec, dDvec; // vectorized edge attributes for solver
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    void regularLattice(const float width, const float height);
    std::vector<unsigned int> rectangularBoundaryIndices();
    std::vector<unsigned int> randomSources(const std::vector<unsigned int>& boundary, unsigned int n);
public:
    Graph(uint32_t seed, const float width, const float height, const unsigned int resolution)
    : m_master_seed { seed }
    , m_rng_sources(seed + 1)
    , m_rng_initD(seed + 2)
    , m_resolution { resolution }
    , I0 { 2.0 * static_cast<double>(resolution) / 4.0 }
    {
        std::cout << "Graph init seed : " << m_master_seed << '\n';
        // vector reservations occur depending on graph initialization type
        regularLattice(width, height);
        initLaplacian();
    };
    
    std::size_t nodeCount() { return m_nodes.size(); }
    std::size_t edgeCount() { return m_edges.size(); }
    unsigned int resolution() const { return m_resolution; }
    const std::vector<Node>& nodes() { return m_nodes; }
    const std::vector<Edge>& edges() { return m_edges; }
    const Eigen::SparseMatrix<double>& getL() { return L; }
    const Eigen::VectorXd& getS() { return s; }
    const Eigen::VectorXd& getD() { return Dvec; }
    const Eigen::VectorXd& getQ() { return Qvec; }
    bool fitConverged() { return fitnessConverged; }
    
    void initLaplacian();
    void updateLaplacian();
    void solvePressures();
    void setSources();
    void computeFlows(bool checkConvergence);
    void updateConductances(double dt);
    double efficiency(const Eigen::VectorXd& D);
    double dissipation(const Eigen::VectorXd& D);
    double transportCost();
    bool conductanceConverged() const;
    bool efficiencyConverged();
    Eigen::VectorXd createPerturbationVec(std::mt19937& rng, double eps);
    double probeHessian(std::mt19937& rng, double eps);
    double probePrune(unsigned int idx, double eps);
    std::vector<double> sampleHSpec(unsigned int nSamples, double eps);

    void solveStep(bool checkConvergence = true)
    {
        updateLaplacian();
        solvePressures(); // automatically sets p
        computeFlows(checkConvergence);
    }

    void evolveGraph(double dt)
    {
        solveStep();
        updateConductances(dt);
    }

    void printSpec(const std::vector<double>& eigvals)
    {
        for (double lambda : eigvals)
        {
            std::cout << lambda << '\n';
        }
    }

};
