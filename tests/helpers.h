#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>

// ElasticRods
#include <ElasticRods/ElasticRod.hh>
#include <ElasticRods/PeriodicRod.hh>
#include <ElasticRods/RodMaterial.hh>

// ElasticKnots
#include "../PeriodicRodList.hh"
#include "../ContactProblem.hh"

PeriodicRod define_periodic_rod(std::vector<Eigen::Vector3d> pts, RodMaterial material);
std::vector<Eigen::Vector3d> read_nodes_from_file(std::string &filename);
