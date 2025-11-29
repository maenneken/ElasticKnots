#include "KnotVisualizer.h"

KnotVisualizer::KnotVisualizer() {
    polyscope::init();
    polyscope::options::autocenterStructures = true;
    polyscope::options::autoscaleStructures = true;
}

void KnotVisualizer::setKnot(const std::vector<Eigen::Vector3d>& pts, const double radius) {
    std::vector<std::array<size_t,2>> edges;
    for (size_t i = 0; i + 1 < pts.size(); ++i)
        edges.push_back({i, i+1});

    curve = polyscope::registerCurveNetwork("knot", pts, edges);
    curve->setRadius(radius, true);
    
}

void KnotVisualizer::updateKnot(const std::vector<Eigen::Vector3d>& pts) {
    if (!curve) return;

    curve->updateNodePositions(pts);
}

void KnotVisualizer::show() {
    polyscope::show();
}
