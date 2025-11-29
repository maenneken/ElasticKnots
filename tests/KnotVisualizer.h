#pragma once
#include <polyscope/curve_network.h>
#include <Eigen/Dense>
#include <vector>

class KnotVisualizer {
public:
    KnotVisualizer();
    void setKnot(const std::vector<Eigen::Vector3d>& pts, const double radius);
    void updateKnot(const std::vector<Eigen::Vector3d>& pts);
    void show();

private:
    polyscope::CurveNetwork* curve = nullptr;
};
