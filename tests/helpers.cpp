#include "helpers.h"
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


/** from helpers.py
    def define_periodic_rod(pts, material, rest_curv_rad=np.inf, total_opening_angle=0, minimize_twist=False):
    duplicated_0 = np.linalg.norm(pts[0, :] - pts[-2, :]) < 1e-12
    duplicated_1 = np.linalg.norm(pts[1, :] - pts[-1, :]) < 1e-12
    if not duplicated_0 and not duplicated_1:
        pts = np.vstack((pts, pts[0, :], pts[1, :]))
    elif duplicated_0 != duplicated_1:
        raise ValueError("Only one of the first two nodes was duplicated.")
        
    pr = elastic_rods.PeriodicRod(pts, zeroRestCurvature=True)  # always set rest curvature to zero, then eventually modify restKappas
    pr.setMaterial(material)
    pr.totalOpeningAngle = total_opening_angle
    
    # Set rest curvature
    if rest_curv_rad != np.inf:
        rest_lengths = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        rest_kappas = compute_rest_kappas(rest_curv_rad=rest_curv_rad, rest_lengths=rest_lengths)
        pr.rod.setRestKappas(rest_kappas)
        
    # Set the bending energy type to match the definition from [Bergou et al. 2010]
    # The bending energy in [Bergou et al. 2008] is technically non-physical.
    pr.rod.bendingEnergyType = elastic_rods.BendingEnergyType.Bergou2010
    
    if minimize_twist:
        # Minimize the twisting energy 
        # (i.e. run an optimization on the \theta variables only, 
        # leaving the ends of the rod free to untwist)
        elastic_knots.minimize_twist(pr)
    
    return pr
**/
PeriodicRod define_periodic_rod(std::vector<Eigen::Vector3d> pts, RodMaterial material){
    std::size_t n_pts = pts.size();
    bool duplicated_0 = (pts[0] - pts[n_pts - 2]).norm()  < 1e-12;
    bool duplicated_1 = (pts[1] - pts[n_pts - 1]).norm()  < 1e-12;
    
    if (!duplicated_0 && !duplicated_1){
        pts.push_back(pts[0]);
        pts.push_back(pts[1]);
    }
    else if (duplicated_0 != duplicated_1){
        throw "Only one of the first two nodes was duplicated.";
    }
    PeriodicRod pr(pts,true);
    pr.setMaterial(material);
    pr.setTotalOpeningAngle(0);
    pr.rod.setBendingEnergyType(ElasticRod::BendingEnergyType::Bergou2010);

    return pr;


}

/**
    def read_nodes_from_file(file):
    """
    Supported extensions: obj, txt
    """
    nodes = []
    connectivity = []
    n_rods = 0
    if file.endswith('.obj'):
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                if line.startswith('v'):
                    pt = []
                    for coord in line.split(' ')[1:4]:
                        pt.append(float(coord))
                    nodes.append(np.array(pt))
                if line.startswith('l'):
                    edge = []
                    for index in line.split(' ')[1:3]:
                        edge.append(int(index))
                        if len(edge) == 2 and int(index) < edge[0]: # last edge of a rod 
                            n_rods += 1
                    connectivity.append(edge)
        if n_rods == 1:
            return np.array(nodes)
        elif n_rods > 1:
            indices_connections = [i for i in range(len(connectivity)) if connectivity[i][0] > connectivity[i][1]]
            ne_per_rod = np.append(indices_connections[0] + 1, np.diff(indices_connections))
            pts = np.array(nodes)
            pts_list = [pts[0:ne_per_rod[0], :]]
            for ri in range(0, n_rods-1):
                pts_list.append(pts[indices_connections[ri]:indices_connections[ri+1]])
            return pts_list
    
    elif file.endswith('.txt'):
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                pt = []
                for coord in line.split(' ')[0:3]:
                    pt.append(float(coord))
                nodes.append(np.array(pt))
        return np.array(nodes)
    
    elif not '.' in file.split('/')[0]: # no extension, assum same formatting as .txt
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                pt = []
                for coord in line.split(' ')[0:3]:
                    pt.append(float(coord))
                nodes.append(np.array(pt))
        return np.array(nodes)
 **/
std::vector<Eigen::Vector3d> read_nodes_from_file(std::string &filename){
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + filename);
    
    std::vector<Eigen::Vector3d> nodes;
    std::string line;

    if (filename.find(".obj") != std::string::npos){ //it is a .obj file
        while (std::getline(file, line)) {
            if (line.empty()) continue;

            std::istringstream stream(line);
            double x, y, z;
            char c;
            stream >> c;
            if(c == 'v' ){ //vertex
                stream >> x >> y >> z;
                nodes.emplace_back(x, y, z);
            }
            else if (c == 'l'){ // connectivity
                continue;   
            }
            else{
                continue; // skip invalid lines
            }
                
        }
        return nodes;
    }
    else if (filename.find(".txt") != std::string::npos){ //it is a .txt file
        throw "not implemented yet";
        return nodes;
    }


}
