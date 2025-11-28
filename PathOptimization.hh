#ifndef SLIDING_PROBLEM_HH
#define SLIDING_PROBLEM_HH

#include <utility>
#include <algorithm>
#include <Eigen/Core>
#include "PeriodicRodList.hh"
#include "SoftConstraint.hh"
#include <tight_inclusion/ccd.hpp>
#include <ipc/ipc.hpp>
#include <MeshFEM/newton_optimizer/newton_optimizer.hh>
#include "ContactProblem.hh"

#include <functional>

using CallbackFunction = std::function<void(NewtonProblem &, size_t)>;

struct PathProblem : public ContactProblem {
    using SoftConstraintsList = std::vector<std::shared_ptr<SoftConstraint>>;

    PathProblem(
        PeriodicRodList &rods, 
        ContactProblemOptions options = ContactProblemOptions()
        ) : m_rods(rods), m_options(options) {

        updateCachedVars();
        updateCachedSparsityPattern();
        updateCharacteristicLength();

        rods.updateSourceFrame();

        if (m_options.hasCollisions) {
            // Initialize CollisionMesh with the connectivity information required by IPC
            Eigen::MatrixXd vertices(m_rods.numVertices(), 3);
            Eigen::MatrixXi edges(m_rods.numEdges(), 2);
            Eigen::MatrixXi faces(0, 3);  // no faces
            for (size_t ri = 0; ri < m_rods.size(); ri++) {
                const size_t fni = m_rods.firstGlobalNodeIndexInRod(ri);
                const size_t nvr = m_rods.numVerticesInRod(ri);
                for (size_t i = 0; i < nvr; i++) {
                    size_t ni = fni + i;
                    vertices.row(ni) = m_rods.getNode(ni);
                    edges(ni, 0) = ni;
                    edges(ni, 1) = ni + 1;
                }
                edges(fni + nvr - 1, 1) = fni;  // overwrite second node in last edge
            }
            m_collisionMesh = ipc::CollisionMesh(vertices, edges, faces);

            updateConstraintSet();
            updateCachedSparsityPattern();
        }

        // Compute minimium edge rest length and check its compatibility with the constraint barrier thickness
        Real min_rl = std::numeric_limits<Real>::infinity();
        for (size_t ri = 0; ri < m_rods.size(); ri++) {
            const std::vector<Real> rli = m_rods[ri]->restLengths();
            Real min_rli = *std::min_element(rli.begin(), rli.end());
            if (min_rli < min_rl)
                min_rl = min_rli;
        }
        if (min_rl < m_options.dHat)
            std::cerr << "WARNING: The minimum edge rest length is smaller than the constraint barrier thickness.\n"
                         "The simulations will continue, but the result might be non-physical due to spurious contact forces between neighboring edges.\n"
                         "Consider decreasing the cross-section radius or using a coarser polyline.\n"
                         "Increasing the minContactEdgeDist parameter would remove non-physical contact forces, too, but only at the expense of topology preservation guarantees.\n" << std::endl;
    }

    virtual const Eigen::VectorXd getVars() const override { return m_cachedVars; }
    virtual size_t numVars() const override { return m_rods.numDoF(); }
    size_t numRods() const { return m_rods.size(); }
    size_t numIPCConstraints() const { return m_constraintSet.size(); }
    
    virtual bool hasCollisions() const { return m_options.hasCollisions; };

    void addSoftConstraint(const std::shared_ptr<SoftConstraint> &sc) { m_softConstraints.push_back(sc); }
    void addSoftConstraints(const std::vector<std::shared_ptr<SoftConstraint>> &scList) { for (auto &sc : scList) m_softConstraints.push_back(sc); }

    //pathEnergy = sum(RodEngery)+sum(distance(rod[i],rod[i+1]))
    virtual Real pathEnergy() const override { 
        Real e = energy();
        Real distance = 0;
        for (size_t i = 0; i < m_rods.size()-1; i++){
            Eigen::VectorXd diff = m_rods[i]->getDoFs() - m_rods[i + 1]->getDoFs();
            distance += diff.squaredNorm();
        }
        return e + distance

    }
    virtual Eigen::VectorXd gradient(bool freshIterate = false) const override {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(numVars());
        result.head(numVars()) = m_rods.gradient(freshIterate);        // rods
        if (m_options.hasCollisions) {
            Eigen::VectorXd bpGrad = m_options.contactStiffness * compute_barrier_potential_gradient(m_collisionMesh, m_rods.deformedPointsMatrix(), m_constraintSet, m_options.dHat);
            for (size_t ri = 0; ri < m_rods.size(); ri++)
                result.segment(m_rods.firstGlobalDofIndexInRod(ri), 3*m_rods.numVerticesInRod(ri)) += bpGrad.segment(3*m_rods.firstGlobalNodeIndexInRod(ri), 3*m_rods.numVerticesInRod(ri));
        }
        if (external_forces.size() > 0) {                              // external potential energy
            assert((size_t)external_forces.size() == numVars());
            result.head(numVars()) -= external_forces;
        }
        for (auto &sc : m_softConstraints)                             // soft constraints
            sc->gradient(m_rods, result);
        return result;
    }
    Real hessianShift = 0.0;  // Multiple of the identity to add to the Hessian on each evaluation.

    mutable SuiteSparseMatrix m_hessianSparsity, m_rodHessianSparsity, m_contactHessianSparsity;
    Real m_characteristicLength = 1.0;
    CallbackFunction m_customCallback;

    PeriodicRodList &m_rods;
    Eigen::VectorXd m_cachedVars;    // [r1, ..., rn], where ri = [x1, y1, z1, ..., xk, yk, zk, th1, ..., thk-1].
    SoftConstraintsList m_softConstraints;
    ContactProblemOptions m_options;

    // ipc-toolkit
    ipc::Constraints m_constraintSet;
    ipc::CollisionMesh m_collisionMesh;  // Mesh connectivity
};


void minimize_twist(PeriodicRod &rod, bool verbose = false);
void spread_twist_preserving_link(PeriodicRod &pr, bool verbose = false);


ConvergenceReport compute_equilibrium(
    PeriodicRodList rods,
    const ContactProblemOptions &problemOptions = ContactProblemOptions(),
    const NewtonOptimizerOptions &optimizerOptions = NewtonOptimizerOptions(), 
    std::vector<size_t> fixedVars = std::vector<size_t>(), 
    const Eigen::VectorXd &externalForces = Eigen::VectorXd(),
    const ContactProblem::SoftConstraintsList &softConstraints = ContactProblem::SoftConstraintsList(),
    CallbackFunction customCallback = nullptr,
    double hessianShift = 0.0
);

#endif /* end of include guard: SLIDING_PROBLEM_HH */