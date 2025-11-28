#include "PathOptimization.hh"

ConvergenceReport compute_equilibrium(
    PeriodicRodList rods,
    const ContactProblemOptions &problemOptions,
    const NewtonOptimizerOptions &optimizerOptions, 
    std::vector<size_t> fixedVars, 
    const Eigen::VectorXd &externalForces,
    const ContactProblem::SoftConstraintsList &softConstraints,
    CallbackFunction customCallback,
    Real hessianShift
    ) {
    std::unique_ptr<ContactProblem> problem = std::make_unique<ContactProblem>(rods, problemOptions);
    problem->addFixedVariables(fixedVars);
    if (externalForces.size() > 0) {
        assert((size_t)externalForces.size() == rods.numDoF());
        problem->external_forces = externalForces;
    }
    problem->addSoftConstraints(softConstraints);
    problem->hessianShift = hessianShift;
    if (customCallback)
        problem->setCustomIterationCallback(customCallback);
    std::unique_ptr<NewtonOptimizer> optimizer = std::make_unique<NewtonOptimizer>(std::move(problem));
    optimizer->options = optimizerOptions;
    return optimizer->optimize();
}