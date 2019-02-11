#ifndef _A_STAR_H
#define _A_STAR_H

#include <ros/console.h>
#include <ros/ros.h>
#include <Eigen/Eigen>
#include <iostream>
// #include "backward.hpp"
#include "grad_traj_optimization/data_type.h"
#include "grad_traj_optimization/sdf_map.h"

class AStarPathFinder {
 private:
  Eigen::Vector3d gridIndex2coord(Eigen::Vector3i index);
  Eigen::Vector3i coord2gridIndex(Eigen::Vector3d pt);
  GridNodePtr pos2gridNodePtr(Eigen::Vector3d pos);

  double getDiagHeu(GridNodePtr node1, GridNodePtr node2);
  double getManhHeu(GridNodePtr node1, GridNodePtr node2);
  double getEuclHeu(GridNodePtr node1, GridNodePtr node2);
  double getHeu(GridNodePtr node1, GridNodePtr node2);

  std::vector<GridNodePtr> retrievePath(GridNodePtr current);

  double resolution, inv_resolution;
  Eigen::Vector3i grid_num;
  Eigen::Vector3d origin;
  double tie_breaker = 1.0 + 1.0 / 10000;

  std::vector<GridNodePtr> expandedNodes;
  std::vector<GridNodePtr> gridPath;

  GridNodePtr*** GridNodeMap;
  std::multimap<double, GridNodePtr> openSet;

 public:
  AStarPathFinder(){};
  ~AStarPathFinder(){};

  void initGridNodeMap(Eigen::Vector3i grid_num, double _resolution,
                       Eigen::Vector3d origin);
  void linkLocalMap(SDFMap& sdf_map);
  bool searchPath(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);

  void resetNode();
  void resetPath();

  std::vector<Eigen::Vector3d> getPath();
  std::vector<GridNodePtr> getVisitedNodes();
};

#endif