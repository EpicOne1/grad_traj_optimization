#include "grad_traj_optimization/a_star.h"

using namespace std;
using namespace Eigen;

void AStarPathFinder::initGridNodeMap(Eigen::Vector3i grid_num,
                                      double _resolution, Vector3d origin) {
  this->origin = origin;
  this->grid_num = grid_num;
  resolution = _resolution;
  inv_resolution = 1.0 / _resolution;

  GridNodeMap = new GridNodePtr**[grid_num(0)];
  for (int i = 0; i < grid_num(0); i++) {
    GridNodeMap[i] = new GridNodePtr*[grid_num(1)];
    for (int j = 0; j < grid_num(1); j++) {
      GridNodeMap[i][j] = new GridNodePtr[grid_num(2)];
      for (int k = 0; k < grid_num(2); k++) {
        Vector3i tmpIdx(i, j, k);
        Vector3d pos = gridIndex2coord(tmpIdx);
        GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
      }
    }
  }
}

void AStarPathFinder::linkLocalMap(SDFMap& sdf_map) {
  Vector3d coord;
  for (int i = 0; i < grid_num(0); i++)
    for (int j = 0; j < grid_num(1); j++) {
      for (int k = 0; k < grid_num(2); k++) {
        GridNodePtr ptr = GridNodeMap[i][j][k];
        ptr->id = 0;
        ptr->occupancy = double(sdf_map.getOccupancy(Eigen::Vector3i(i, j, k)));
        ptr->distance = double(sdf_map.getDistance(Eigen::Vector3i(i, j, k)));
      }
    }
}

void AStarPathFinder::resetNode() {
  // ROS_WARN("expandedNodes size : %d", expandedNodes.size());
  for (auto tmpPtr : expandedNodes) {
    tmpPtr->occupancy = 0;  // forget the occupancy
    tmpPtr->id = 0;
    tmpPtr->cameFrom = NULL;
    tmpPtr->gScore = inf;
    tmpPtr->fScore = inf;
  }

  for (auto ptr : openSet) {
    GridNodePtr tmpPtr = ptr.second;
    tmpPtr->occupancy = 0;  // forget the occupancy
    tmpPtr->id = 0;
    tmpPtr->cameFrom = NULL;
    tmpPtr->gScore = inf;
    tmpPtr->fScore = inf;
  }

  expandedNodes.clear();
  // ROS_WARN("local map reset finish");
}

GridNodePtr AStarPathFinder::pos2gridNodePtr(Vector3d pos) {
  Vector3i idx = coord2gridIndex(pos);
  GridNodePtr grid_ptr = new GridNode(idx, pos);

  return grid_ptr;
}

Vector3d AStarPathFinder::gridIndex2coord(Vector3i index) {
  Vector3d pt;
  // cell_x_size_ * ((double)x_index + 0.5), cell_y_size_ * ((double)y_index +
  // 0.5), cell_z_size_ * ((double)z_index + 0.5)

  pt(0) = ((double)index(0) + 0.5) * resolution + origin(0);
  pt(1) = ((double)index(1) + 0.5) * resolution + origin(1);
  pt(2) = ((double)index(2) + 0.5) * resolution + origin(2);

  /*pt(0) = (double)index(0) * resolution + origin(0) + 0.5 * resolution;
  pt(1) = (double)index(1) * resolution + origin(1) + 0.5 * resolution;
  pt(2) = (double)index(2) * resolution + origin(2) + 0.5 * resolution;*/
  return pt;
}

Vector3i AStarPathFinder::coord2gridIndex(Vector3d pt) {
  Vector3i idx;
  idx << min(max(int((pt(0) - origin(0)) * inv_resolution), 0),
             grid_num(0) - 1),
      min(max(int((pt(1) - origin(1)) * inv_resolution), 0), grid_num(1) - 1),
      min(max(int((pt(2) - origin(2)) * inv_resolution), 0), grid_num(2) - 1);

  return idx;
}

double AStarPathFinder::getDiagHeu(GridNodePtr node1, GridNodePtr node2) {
  double dx = abs(node1->index(0) - node2->index(0));
  double dy = abs(node1->index(1) - node2->index(1));
  double dz = abs(node1->index(2) - node2->index(2));

  double h;
  int diag = min(min(dx, dy), dz);
  dx -= diag;
  dy -= diag;
  dz -= diag;

  if (dx == 0) {
    h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dy, dz) + 1.0 * abs(dy - dz);
  }
  if (dy == 0) {
    h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dz) + 1.0 * abs(dx - dz);
  }
  if (dz == 0) {
    h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dy) + 1.0 * abs(dx - dy);
  }
  return h;
}

double AStarPathFinder::getManhHeu(GridNodePtr node1, GridNodePtr node2) {
  double dx = abs(node1->index(0) - node2->index(0));
  double dy = abs(node1->index(1) - node2->index(1));
  double dz = abs(node1->index(2) - node2->index(2));

  return dx + dy + dz;
}

double AStarPathFinder::getEuclHeu(GridNodePtr node1, GridNodePtr node2) {
  return (node2->index - node1->index).norm();
}

double AStarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2) {
  return tie_breaker * getDiagHeu(node1, node2);
  // return tie_breaker * getEuclHeu(node1, node2);
}

vector<GridNodePtr> AStarPathFinder::retrievePath(GridNodePtr current) {
  vector<GridNodePtr> path;
  path.push_back(current);

  while (current->cameFrom != NULL) {
    current = current->cameFrom;
    path.push_back(current);
  }

  return path;
}

vector<GridNodePtr> AStarPathFinder::getVisitedNodes() {
  vector<GridNodePtr> d;
  for (int i = 0; i < grid_num(0); i++)
    for (int j = 0; j < grid_num(1); j++)
      for (int k = 0; k < grid_num(2); k++) {
        if (GridNodeMap[i][j][k]->id != 0)
          // if(GridNodeMap[i][j][k]->id == -1)
          d.push_back(GridNodeMap[i][j][k]);
      }

  cout << "visited node size:" << d.size() << endl;
  return d;
}

/*bool AStarPathFinder::minClearance()
{
    neighborPtr->occupancy > 0.5
}
*/
bool AStarPathFinder::searchPath(Eigen::Vector3d start_pt,
                                 Eigen::Vector3d end_pt) {
  ros::Time time_1 = ros::Time::now();

  Vector3i start_idx = coord2gridIndex(start_pt);
  Vector3i end_idx = coord2gridIndex(end_pt);
  GridNodePtr startPtr = GridNodeMap[start_idx(0)][start_idx(1)][start_idx(2)];
  GridNodePtr endPtr = GridNodeMap[end_idx(0)][end_idx(1)][end_idx(2)];

  openSet.clear();

  GridNodePtr neighborPtr = NULL;
  GridNodePtr current = NULL;

  startPtr->gScore = 0;
  startPtr->fScore = getHeu(startPtr, endPtr);
  startPtr->id = 1;  // put start node in open set
  startPtr->coord = start_pt;
  openSet.insert(
      make_pair(startPtr->fScore, startPtr));  // put start in open set

  double tentative_gScore;

  int num_iter = 0;
  while (!openSet.empty()) {
    num_iter++;
    current = openSet.begin()->second;

    if (current->index(0) == endPtr->index(0) &&
        current->index(1) == endPtr->index(1) &&
        current->index(2) == endPtr->index(2)) {
      ROS_WARN("[Astar]Reach goal..");
      // cout << "goal coord: " << endl << current->real_coord << endl;
      cout << "total number of iteration used in Astar: " << num_iter << endl;
      ros::Time time_2 = ros::Time::now();
      ROS_WARN("Time consume in A star path finding is %f",
               (time_2 - time_1).toSec());
      gridPath = retrievePath(current);
      return true;
    }
    openSet.erase(openSet.begin());
    current->id = -1;  // move current node from open set to closed set.
    expandedNodes.push_back(current);

    for (int dx = -1; dx < 2; dx++)
      for (int dy = -1; dy < 2; dy++)
        for (int dz = -1; dz < 2; dz++) {
          if (dx == 0 && dy == 0 && dz == 0) continue;

          Vector3i neighborIdx;
          neighborIdx(0) = (current->index)(0) + dx;
          neighborIdx(1) = (current->index)(1) + dy;
          neighborIdx(2) = (current->index)(2) + dz;

          if (neighborIdx(0) < 0 || neighborIdx(0) >= grid_num(0) ||
              neighborIdx(1) < 0 || neighborIdx(1) >= grid_num(1) ||
              neighborIdx(2) < 0 || neighborIdx(2) >= grid_num(2)) {
            continue;
          }

          neighborPtr =
              GridNodeMap[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];

          /*                    if(minClearance() == false){
                                  continue;
                              }*/

          if (neighborPtr->distance < 0.5) {
            continue;
          }

          if (neighborPtr->id == -1) {
            continue;  // in closed set.
          }

          double static_cost = sqrt(dx * dx + dy * dy + dz * dz);

          tentative_gScore = current->gScore + static_cost;

          if (neighborPtr->id != 1) {
            // discover a new node
            neighborPtr->id = 1;
            neighborPtr->cameFrom = current;
            neighborPtr->gScore = tentative_gScore;
            neighborPtr->fScore =
                neighborPtr->gScore + getHeu(neighborPtr, endPtr);
            neighborPtr->nodeMapIt = openSet.insert(make_pair(
                neighborPtr->fScore,
                neighborPtr));  // put neighbor in open set and record it.
            continue;
          } else if (tentative_gScore <=
                     neighborPtr->gScore) {  // in open set and need update
            neighborPtr->cameFrom = current;
            neighborPtr->gScore = tentative_gScore;
            neighborPtr->fScore =
                tentative_gScore + getHeu(neighborPtr, endPtr);
            openSet.erase(neighborPtr->nodeMapIt);
            neighborPtr->nodeMapIt = openSet.insert(make_pair(
                neighborPtr->fScore,
                neighborPtr));  // put neighbor in open set and record it.
          }
        }
  }

  ros::Time time_2 = ros::Time::now();
  ROS_WARN("Time consume in A star path finding is %f",
           (time_2 - time_1).toSec());
  return false;
}

vector<Vector3d> AStarPathFinder::getPath() {
  vector<Vector3d> path;

  for (auto ptr : gridPath) path.push_back(ptr->coord);

  reverse(path.begin(), path.end());
  return path;
}

void AStarPathFinder::resetPath() { gridPath.clear(); }