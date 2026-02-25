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

    // std::vector<unsigned int> boundary { rectangularBoundaryIndices() };
    // m_source_ids = randomSources(boundary, n_sources);

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
    std::normal_distribution<double> noise(0.0, 2e-1); // stdev ~ 0.01% relative perturbation

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
                // adding noise here is OK since |noise| < 1.0 always (no negative D ever)
                // std::cout << (1.0 + noise(m_rng_initD)) << '\n';
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
    // assert(m_sink_idx < s.size());
    setSources();

    sr.resize(N-1);

    Eigen::VectorXi nnzPerRow { Eigen::VectorXi::Constant(N, 5) };
    L.reserve(nnzPerRow);

    updateLaplacian();
    solvePressures();
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

void Graph::solvePressures()
{
    int k { 0 }; // pick a zero-pressure node
    // int k { static_cast<int>(m_sink_idx) };
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
    // solver.analyzePattern(Lr);
    solver.factorize(Lr);

    if (solver.info() == Eigen::Success)
    {
        Eigen::VectorXd pr = solver.solve(sr);
        // reconstruct full p
        p.setZero();
        for (int i = 0; i < N; ++i)
            if (i != k)
                p(i) = pr(map(i));
    }
    else
    {
        // Handle the error (e.g., matrix is singular)
        std::cerr << "Decomposition Failed: " << solver.lastErrorMessage() << std::endl;
    }
}

void Graph::setSources()
{
    s.setZero();
    // s(m_sink_idx) = -I0 * static_cast<double>(m_source_ids.size());
    // // double Isrc { I0 / static_cast<double>(m_source_ids.size()) };
    // for (unsigned int idx : m_source_ids) { s(idx) = I0; }
    
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_int_distribution<unsigned int> dist(0, static_cast<unsigned int>(s.size())-1);

    // for (unsigned int i = 0; i < s.size(); ++i)
    for (unsigned int i = 0; i < fmin(n_sources, s.size()); ++i)
    {
        // s(i) = normal(m_rng_sources);
        double s_i { abs(normal(m_rng_sources)) };
        unsigned int idx { dist(m_rng_sources) };
        if (idx != m_sink_idx)
        {
            s(dist(m_rng_sources)) += s_i;
        }
        // s(dist(m_rng_sources)) += s_i;
        // s(dist(m_rng_sources)) -= s_i;
    }

    // s.array() -= s.mean();
    s(m_sink_idx) = -s.sum();
    s.normalize();
    s *= I0;
}

void Graph::computeFlows(bool checkConvergence)
{
    // std::cout << "###############################" << '\n';
    // for dissipation energy metric
    double E_old { dissipation(Dvec-dDvec) };
    // double E_old { E };
    // E = 0.0;

    for (unsigned int k = 0; k < Dvec.size(); ++k)
    {
        const Edge& edge { m_edges[k] };
        
        Qvec(k) = Dvec(k) * (p(edge.i) - p(edge.j));
        // std::cout << "D_k : " << '\t' << Dvec(k) << '\t' << "|Q_k| : " << '\t' << abs(Qvec(k)) << '\t' << "E_k : " << '\t' << Qvec(k) * Qvec(k) / Dvec(k) << '\n';
        // E += Qvec(k) * Qvec(k) / Dvec(k); // might need its own wrapper but fine for now
    }

    // check fitness (dissipation) for convergence
    if (checkConvergence && ((dissipation(Dvec) - E_old) / E_old < m_tol)) { fitnessConverged = true; }
    // if (checkConvergence && (abs(E - E_old) / E_old < m_tol)) { fitnessConverged = true; }
}

void Graph::updateConductances(const double dt)
{
    constexpr double eps { 1e-12 }; // 1e-16 too small --> leads to singular L during perturbations
    double alpha { 100.0 };
    double beta { 10.0 };
    double gamma { 3.0 };

    for (unsigned int k = 0; k < Dvec.size(); ++k)
    {
        double Qgamma { pow(abs(Qvec(k)), gamma) };
        dDvec(k) = dt * ( alpha * Qgamma / (1.0 + Qgamma) - beta * Dvec(k) );
        // dDvec(k) = dt * ( alpha * pow(abs(Qvec(k)), gamma) - beta * Dvec(k) );
        Dvec(k) += dDvec(k);
        Dvec(k) = fmax(Dvec(k), eps);
    }
}

double Graph::dissipation(const Eigen::VectorXd& D)
{
    // assumes conductances have already been updated to next step
    return (Qvec.cwiseProduct(Qvec).cwiseQuotient(D) + c_t * D.cwisePow(0.5)).sum();
}

double Graph::efficiency(const Eigen::VectorXd& D)
{
    return s(m_sink_idx) * s(m_sink_idx) / D.sum();
}

double Graph::transportCost()
{
    return s.transpose() * L.adjoint() * s;
}

bool Graph::conductanceConverged() const
{
    return dDvec.norm() / Dvec.norm() < m_tol;
}

Eigen::VectorXd Graph::createPerturbationVec(std::mt19937& rng, double eps, std::vector<unsigned int> aliveIdxs)
{
    std::normal_distribution<double> noise(0.0, 1.0);
    // std::normal_distribution<double> noise(1.0, eps);
    unsigned int N { static_cast<unsigned int>(Dvec.size()) };
    Eigen::VectorXd v(N);
    // Eigen::VectorXd v { Eigen::VectorXd::Zero(N) };

    for (unsigned int i = 0; i < N; ++i)
    // for (unsigned int i : aliveIdxs)
    {
        // Rademacher generator
        // v(i) = (rng() & 1) ? 1 : -1;
        // v(i) = (rng() & 1) ? exp(eps * noise(rng)) : exp(-eps * noise(rng));
        // v(i) = exp(eps * noise(rng) / sqrt(N)); // norm approximation is roughly okay, especially for large N
        
        v(i) = noise(rng);
    }

    // return v / sqrt(N);
    // return eps * v.normalized();
    return (eps * v.normalized()).array().exp();
}

double Graph::probeHessian(std::mt19937& rng, double eps, std::vector<unsigned int> aliveIdxs)
{
    // store converged fitness
    // double Fstar { efficiency(Dvec) };
    // double Fstar { E };
    solveStep(false);
    double Fstar { dissipation(Dvec) };
    
    Eigen::VectorXd delta { createPerturbationVec(rng, eps, aliveIdxs) };
    // delta.fill(eps); // to perturb in Dvec-direction explicitly
    // for (int i = 0; i < delta.size(); ++i)
    // {
    //     std::cout << "delta_" << i << " = " << '\t' << delta(i) << '\n';
    // }
    // Eigen::VectorXd one(delta.size());
    // one.fill(1.0);

    // since ||D|| grows with network size (roughly sqrt(dim(D)) * avg D_i ),
    // perturbation is scaled by the dimensionless parameter ||D||/sqrt(dim(D))
    // std::cout << "|delta| = " << eps * Dvec.norm() / sqrt(Dvec.size()) << '\n';
    // delta *= eps * Dvec.norm() / sqrt(Dvec.size());

    // copy current state
    Eigen::VectorXd Dstar = Dvec;
    Eigen::VectorXd Qstar = Qvec;
    
    // finite difference calculation requires Dvec^* +/- delta
    // + delta
    // Dvec += delta;
    // Eigen::VectorXd opd { one + delta };
    // for (int i = 0; i < delta.size(); ++i)
    // {
    //     std::cout << "D_" << i << " = " << '\t' << Dstar(i) << '\n';
    //     std::cout << "(one + delta)_" << i << " = " << '\t' << opd(i) << '\n';
    // }
    // Dvec = Dstar.cwiseProduct(opd);
    Dvec = Dstar.cwiseProduct(delta);
    
    // std::cout << "Relative D+ : " << '\t' << Dvec.norm() / Dstar.norm() << '\n';
    
    // recompute flows to get fitness
    solveStep(false);
    // double Fplus { efficiency(Dvec) };
    // double Fplus { E };
    double Fplus { dissipation(Dvec) };
    
    // - delta
    // Dvec -= 2 * delta;
    // Eigen::VectorXd omd { one - delta };
    // for (int i = 0; i < delta.size(); ++i)
    // {
    //     std::cout << "D_" << i << " = " << '\t' << Dstar(i) << '\n';
    //     std::cout << "(one - delta)_" << i << " = " << '\t' << omd(i) << '\n';
    // }
    // Dvec = Dstar.cwiseProduct(omd);
    Dvec = Dstar.cwiseQuotient(delta);
    // std::cout << "Relative D- : " << '\t' << Dvec.norm() / Dstar.norm() << '\n';
    solveStep(false);
    // double Fminus { efficiency(Dvec) };
    // double Fminus { E };
    double Fminus { dissipation(Dvec) };

    // revert system to steady state
    // (don't need to recompute flows unless the simulation will continue to evolve beyond Hessian probing)
    Dvec = Dstar - dDvec;
    Qvec = Qstar;
    solveStep(false);
    Dvec = Dstar;

    // return (Fplus - 2.0 * Fstar + Fminus) / delta.squaredNorm();
    // exp(epsilon) ~ 1 + epsilon ==> D exp(epsilon) ~ D + epsilon D
    // D_i * delta_i ~ D (1 +/- epsilon)
    
    // normalize by Fstar to compile Rayleigh coefficients across multiple graph instances
    if (solver.info() != Eigen::Success)
    {
        return -1;
    }
    else
    {
        return (Fplus - 2.0 * Fstar + Fminus) / (eps * eps) / Fstar; // step size division is ambiguous if use delta from Gaussian
    }
}

std::vector<double> Graph::sampleHSpec([[maybe_unused]] unsigned int nSamples, double eps)
{
    std::mt19937 rng(m_master_seed);
    
    // for probing via edge pruning

    std::vector<unsigned int> aliveEdgeIdx; // (static_cast<unsigned int>(Dvec.size()));
    for (unsigned int i = 0; i < Dvec.size(); ++i)
    {
        if ( (Dvec(i) < 1e-4))
        {
            continue;
        }
        aliveEdgeIdx.push_back(i);
    }

    std::vector<double> eigvals;
    
    // for pruning
    // eigvals.reserve(static_cast<unsigned int>(aliveEdgeIdx.size()));
    // std::cout << aliveEdgeIdx.size() << " edges alive" << '\n';
    // for (unsigned int idx : aliveEdgeIdx)
    // {
    //     eigvals.emplace_back(probePrune(idx, eps));
    // }
    
    // for scaling perturbation
    eigvals.reserve(nSamples);
    for (unsigned int i = 0; i < nSamples; ++i)
    {
        eigvals.push_back(probeHessian(rng, eps, aliveEdgeIdx));
    }

    return eigvals;
}

bool Graph::efficiencyConverged()
{
    double F_old { efficiency(Dvec - dDvec) };
    double F { efficiency(Dvec) };

    return (abs(F - F_old) / F_old) < m_tol;
}

double Graph::probePrune(unsigned int idx, double eps)
{
    solveStep(false);
    double Fstar { dissipation(Dvec) };

    // copy current state
    Eigen::VectorXd Dstar = Dvec;
    Eigen::VectorXd Qstar = Qvec;

    // prune/block edge
    Dvec(idx) *= eps;
    solveStep(false);
    double F { dissipation(Dvec) };

    // revert system to steady state
    Dvec = Dstar - dDvec;
    Qvec = Qstar;
    solveStep(false);
    Dvec = Dstar;

    return (F - Fstar) / Fstar;
}

