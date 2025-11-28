#include "helpers.h"
#include <iostream>
#include <vector>

// ElasticRods
#include <ElasticRods/ElasticRod.hh>
#include <ElasticRods/PeriodicRod.hh>
#include <ElasticRods/RodMaterial.hh>

// ElasticKnots
#include "../PeriodicRodList.hh"
#include "../ContactProblem.hh"


int main() {

    //define Contact Problem
    std::string file = "../data/L400-r0.2-UpTo9Crossings/4_1/0033.obj";
    double rod_radius = 0.2;
    std::vector<double> params = {rod_radius, rod_radius};
    RodMaterial material(
        "ellipse",
        2000,     // Young's modulus
        0.3,      // Poisson's ratio
        params
    );
    std::vector<Eigen::Vector3d> centerline = read_nodes_from_file(file);
    PeriodicRod pr = define_periodic_rod(centerline,material);
    PeriodicRodList rod_list = PeriodicRodList(pr);
    ContactProblemOptions problemOptions;
    problemOptions.hasCollisions = true;
    problemOptions.contactStiffness = 10000;
    problemOptions.dHat = 2 * rod_radius;
    ContactProblem cp(rod_list, problemOptions);

    TripletMatrix<Triplet<Real>> H_sparse = cp.hessian();

    
    Eigen::MatrixXd H;  // symmetric matrix

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);

    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();
    
    return 0;
}
