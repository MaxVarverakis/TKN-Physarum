#include "Graph.hpp"

std::vector<unsigned int> Graph::rectangularBoundaryIndices()
{
    unsigned int N { m_resolution };
    std::vector<unsigned int> indices;
    // 2N + 2(N-2) = 4N - 4 = 4(N-1)
    indices.reserve(4 * (N - 1));

    // left edge
    for (unsigned int col = 0; col < N; ++col)
    {
        indices.push_back(col);
    }
    
    // right edge
    for (unsigned int col = 0; col < N; ++col)
    {
        indices.push_back((N-1) * N + col);
    }

    // bottom edge
    for (unsigned int row = 1; row < N-1; ++row)
    {
        indices.push_back(row * N);
    }
    
    // top edge
    for (unsigned int row = 1; row < N-1; ++row)
    {
        indices.push_back(row * N - 1);
    }

    return indices;
}

std::vector<unsigned int> Graph::randomSources(const std::vector<unsigned int>& boundary, unsigned int n)
{
    if (n > boundary.size()) { n = static_cast<unsigned int>(boundary.size()); }

    std::vector<unsigned int> copy = boundary;

    std::shuffle(copy.begin(), copy.end(), m_rng_sources);
    copy.resize(n); // keep only first n shuffled elements
    return copy;
}

void Graph::regularLattice(const float width, const float height)
{
    m_sink_idx = (m_resolution + 1) * (m_resolution - 1) / 2;

    std::vector<unsigned int> boundary { rectangularBoundaryIndices() };
    m_source_ids = randomSources(boundary, 6);

    // Determine the aspect ratio (with padding)
    float pad { 0.05f };
    float pW { width * (1.0f - 2.0f * pad) };
    float pH { height * (1.0f - 2.0f * pad) };
    
    float dx { pW / (m_resolution - 1) };
    float dy { pH / (m_resolution - 1) };

    m_nodes.reserve( m_resolution * m_resolution );
    std::size_t ec { 2 * m_resolution * (m_resolution - 1) };
    m_edges.reserve(ec); // 2n(n-1)

    // TODO: verify Hessian spectrum histogram is preserved with: (1) no noise, (2) 1e-4 noise, (3) hopefully ok @ 1e-3 noise too
    std::normal_distribution<double> noise(0.0, 1e-4); // stdev ~ 0.01% relative perturbation

    // prepare vectorized versions of edge attributes
    Dvec.resize(static_cast<int>(ec));
    Qvec.resize(static_cast<int>(ec));
    Qvec.setZero();
    dDvec.resize(static_cast<int>(ec));
    dDvec.setZero();

    for (unsigned int col_idx = 0; col_idx < m_resolution; ++col_idx)
    {
        for (unsigned int row_idx = 0; row_idx < m_resolution; ++row_idx)
        {
            // manually add source/sink terms after, otherwise we will need to check if on a source/sink index for every node initialization!
            m_nodes.emplace_back(Node{glm::fvec2(col_idx * dx + width * pad, row_idx * dy + height * pad)});
            
            if (col_idx != m_resolution - 1)
            {
                // right (horizontal) edge
                m_edges.emplace_back(col_idx * m_resolution + row_idx, (col_idx + 1) * m_resolution + row_idx, D0 * (1.0 + noise(m_rng_initD)), 0.0);
            }
            if (row_idx != m_resolution - 1)
            {
                // up (vertical) edge
                m_edges.emplace_back(col_idx * m_resolution + row_idx, col_idx * m_resolution + row_idx + 1, D0 * (1.0 + noise(m_rng_initD)), 0.0);
            }
        }
    }

    // initialize Dvec
    for (unsigned int i = 0; i < ec; ++i)
    {
        Dvec(i) = m_edges[i].D;
    }
};

void Graph::initLaplacian()
{
    int N { static_cast<int>(m_nodes.size()) };
    L.resize(N, N);
    Lr.resize(N-1, N-1);
    
    p.resize(N);

    s.resize(N);
    // merely for visualizing nodes on first frame:
    assert(m_sink_idx < s.size());
    solvePressures();

    sr.resize(N-1);

    Eigen::VectorXi nnzPerRow { Eigen::VectorXi::Constant(N, 5) };
    L.reserve(nnzPerRow);

    updateLaplacian();
    reduceSystem();
}

void Graph::updateLaplacian()
{
    L.setZero();

    for (unsigned int k = 0; k < Dvec.size(); ++k)
    {
        const Edge& edge { m_edges[k] };

        unsigned int i { edge.i };
        unsigned int j { edge.j };
        double D { Dvec(k) };

        L.coeffRef(i, j) -= D;
        L.coeffRef(j, i) -= D;

        L.coeffRef(i, i) += D;
        L.coeffRef(j, j) += D;
    }

    L.makeCompressed();
}

void Graph::reduceSystem()
{
    int k { static_cast<int>(m_sink_idx) };
    int N { static_cast<int>(m_nodes.size()) };

    // index mapping
    auto map = [&](int i){ return (i < k) ? i : i-1; };

    // fill reduced Laplacian
    std::vector<Eigen::Triplet<double>> trips;
    trips.clear();
    trips.reserve(static_cast<std::size_t>(L.nonZeros()));
    for (int n = 0; n < N; ++n)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, n); it; ++it)
        {
            int i = static_cast<int>(it.row());
            int j = static_cast<int>(it.col());
            if (i == k || j == k) continue;
            trips.emplace_back(map(i), map(j), it.value());
        }
    }
    Lr.setFromTriplets(trips.begin(), trips.end());
    Lr.makeCompressed();

    // fill reduced s
    for (int i = 0; i < N; ++i)
    {
        if (i == k) continue;
        sr(map(i)) = s(i);
    }

    if (!solverInitialized) {
        solver.analyzePattern(Lr);
        solverInitialized = true;
    }

    // solve pressures
    solver.factorize(Lr);
    Eigen::VectorXd pr = solver.solve(sr);

    // reconstruct full p
    p.setZero();
    for (int i = 0; i < N; ++i)
        if (i != k)
            p(i) = pr(map(i));
}

void Graph::solvePressures()
{
    s.setZero();
    s(m_sink_idx) = -I0 * static_cast<double>(m_source_ids.size());
    // double Isrc { I0 / static_cast<double>(m_source_ids.size()) };
    for (unsigned int idx : m_source_ids) { s(idx) = I0; }
}

void Graph::computeFlows()
{   
    // for fitness metrics
    double E_old { E };
    totalD = 0.0;
    E = 0.0;

    for (unsigned int k = 0; k < Dvec.size(); ++k)
    {
        const Edge& edge { m_edges[k] };
        
        Qvec(k) = Dvec(k) * (p(edge.i) - p(edge.j));

        totalD += Dvec(k);
        E += Qvec(k) * Qvec(k) / Dvec(k); // might need its own wrapper but fine for now
    }

    // check fitness (dissipation) for convergence
    if (abs(E - E_old) / E_old < m_tol) { fitnessConverged = true; }
}

void Graph::updateConductances(double dt)
{
    constexpr double eps { 1e-8 };
    double alpha { 10.0 };
    double gamma { 2.0 }; // != 1 for nonlinear behavior

    for (unsigned int k = 0; k < Dvec.size(); ++k)
    {
        dDvec(k) = pow(alpha * abs(Qvec(k)), gamma) - Dvec(k);
        Dvec(k) += dt * dDvec(k);
        Dvec(k) = fmax(Dvec(k), eps);
    }
}

double Graph::efficiency()
{
    return s(m_sink_idx) * s(m_sink_idx) / totalD;
}

bool Graph::conductanceConverged() const
{
    return dDvec.norm() / Dvec.norm() < m_tol;
}
