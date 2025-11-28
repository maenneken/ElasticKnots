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
    std::string file = "../data/L400-r0.2-UpTo9Crossings/4_1/0033.obj";
    double rod_radius = 0.2;
    std::vector<double> params = {rod_radius, rod_radius};
    // 1. Create material
    RodMaterial material(
        "ellipse",
        2000,     // Young's modulus
        0.3,      // Poisson's ratio
        params
    );

    // 2. Read centerline nodes
    /** 
    std::vector<Eigen::Vector3d> centerline = {
    {0,0,0},
    {0,1,0},
    {1,1,0},
    {1,0,0},
    {0,0,0},
    {0,1,0}
    };
    **/
    std::vector<Eigen::Vector3d> centerline = read_nodes_from_file(file);

    PeriodicRod pr = define_periodic_rod(centerline,material);

    std::cout << "created pr"<< std::endl;

    // 4. Wrap into PeriodicRodList
    PeriodicRodList rod_list = PeriodicRodList(pr);

    std::cout << "created rodlist"<< std::endl;

    // 5. Setup problem options
    ContactProblemOptions problemOptions;
    problemOptions.hasCollisions = true;
    problemOptions.contactStiffness = 10000;
    problemOptions.dHat = 2 * rod_radius;

    std::cout << "set pr options"<< std::endl;

    // 6. Create contact problem
    ContactProblem cP(rod_list, problemOptions);
    std::cout << "created CP"<< std::endl;

    // 7. Output gradient and DoF sizes
    std::cout << "Gradient length: " << cP.gradient().size()
              << ", DoFs length: " << pr.getDoFs().size()
              << ", same as Vars: " << cP.getVars().size()
              << ", Energy: " << cP.energy()
              << ", contact Energy: " << cP.contactEnergy()
              << std::endl;

    // rod_list.hessian()
    cP.hessian();
    compute_equilibrium(rod_list,problemOptions);
    
    return 0;
}
