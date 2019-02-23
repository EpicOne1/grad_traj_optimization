
#include <grad_traj_optimization/hybrid_astar.h>
#include <sstream>

using namespace std;
using namespace Eigen;

void HybridAStarPathFinder::setParameterAuto() {
  // ros::param::get("/kgb_traj/max_tau", max_tau);
  // ros::param::get("/kgb_traj/init_max_tau", init_max_tau);
  // ros::param::get("/kgb_traj/max_vel", max_vel);
  // ros::param::get("/kgb_traj/max_acc", max_acc);
  // ros::param::get("/kgb_traj/max_jerk", max_jerk);
  // ros::param::get("/kgb_traj/w_time", w_time);
  // ros::param::get("/kgb_traj/horizon", horizon);
  // ros::param::get("/kgb_traj/lambda_heu", lambda_heu);

  max_tau = 0.4;
  max_vel = 3.0;
  max_acc = 2.0;
  w_time = 15.0;
  horizon = 50.0;
  lambda_heu = 5.0;
}


void HybridAStarPathFinder::initGridNodeMap(Eigen::Vector3i grid_num,
                                            double _resolution,
                                            Vector3d origin) {
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

  /* state transition */
  phi.resize(6, 6);
  phi.block(3, 0, 3, 3) = Eigen::MatrixXd::Zero(3, 3);
  phi.block(0, 0, 3, 3) = phi.block(3, 3, 3, 3) =
      Eigen::MatrixXd::Identity(3, 3);
  phi.block(0, 3, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
}

void HybridAStarPathFinder::linkLocalMap(SDFMap& sdf_map) {
  Vector3d coord;
  for (int i = 0; i < grid_num(0); i++)
    for (int j = 0; j < grid_num(1); j++) {
      for (int k = 0; k < grid_num(2); k++) {
        GridNodePtr ptr = GridNodeMap[i][j][k];
        ptr->id = 0;
        ptr->occupancy = double(sdf_map.getOccupancy(Eigen::Vector3i(i, j, k)));
        ptr->distance = sdf_map.getDistance(Eigen::Vector3i(i, j, k));
      }
    }
  // cout << endl;
}

void HybridAStarPathFinder::resetNode() {
  for (auto tmpPtr : expandedNodes) {
    tmpPtr->id = 0;
    tmpPtr->cameFrom = NULL;
    tmpPtr->gScore = inf;
    tmpPtr->fScore = inf;
  }

  for (auto ptr : openSet) {
    GridNodePtr tmpPtr = ptr.second;
    tmpPtr->id = 0;
    tmpPtr->cameFrom = NULL;
    tmpPtr->gScore = inf;
    tmpPtr->fScore = inf;
  }

  expandedNodes.clear();
}

void HybridAStarPathFinder::resetLocalMap() {
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

GridNodePtr HybridAStarPathFinder::pos2gridNodePtr(Vector3d pos) {
  Vector3i idx = coord2gridIndex(pos);
  GridNodePtr grid_ptr = new GridNode(idx, pos);

  Eigen::VectorXd state(6);
  state.head(3) = pos;
  state.tail(3) = Eigen::Vector3d::Zero();
  grid_ptr->state = state;

  return grid_ptr;
}

Vector3d HybridAStarPathFinder::gridIndex2coord(Vector3i index) {
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

Vector3i HybridAStarPathFinder::coord2gridIndex(Vector3d pt) {
  Vector3i idx;
  idx << min(max(int((pt(0) - origin(0)) * inv_resolution), 0),
             grid_num(0) - 1),
      min(max(int((pt(1) - origin(1)) * inv_resolution), 0), grid_num(1) - 1),
      min(max(int((pt(2) - origin(2)) * inv_resolution), 0), grid_num(2) - 1);

  return idx;
}

double HybridAStarPathFinder::getDiagHeu(GridNodePtr node1, GridNodePtr node2) {
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

double HybridAStarPathFinder::getManhHeu(GridNodePtr node1, GridNodePtr node2) {
  double dx = abs(node1->index(0) - node2->index(0));
  double dy = abs(node1->index(1) - node2->index(1));
  double dz = abs(node1->index(2) - node2->index(2));

  return dx + dy + dz;
}

double HybridAStarPathFinder::getEuclHeu(GridNodePtr node1, GridNodePtr node2) {
  return (node2->index - node1->index).norm();
}

double HybridAStarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2) {
  return tie_breaker * getDiagHeu(node1, node2);
  // return tie_breaker * getEuclHeu(node1, node2);
}

vector<GridNodePtr> HybridAStarPathFinder::retrievePath(GridNodePtr current) {
  vector<GridNodePtr> path;
  path.push_back(current);

  while (current->cameFrom != NULL) {
    current = current->cameFrom;
    path.push_back(current);
  }

  return path;
}

vector<GridNodePtr> HybridAStarPathFinder::getVisitedNodes() {
  vector<GridNodePtr> d;
  for (int i = 0; i < grid_num(0); i++)
    for (int j = 0; j < grid_num(1); j++)
      for (int k = 0; k < grid_num(2); k++) {
        if (GridNodeMap[i][j][k]->id != 0)
          // if(GridNodeMap[i][j][k]->id == -1)
          d.push_back(GridNodeMap[i][j][k]);
      }

  // cout << "visited node size:" << d.size() << endl;
  return d;
}

/*bool HybridAStarPathFinder::minClearance()
{
    neighborPtr->occupancy > 0.5
}
*/
int HybridAStarPathFinder::searchPath(Eigen::Vector3d start_pt,
                                      Eigen::Vector3d start_v,
                                      Eigen::Vector3d start_a,
                                      Eigen::Vector3d end_pt,
                                      Eigen::Vector3d end_v, bool init) {
  // cout << "[Hybrid A Star]: begin-----------" << endl;
  // cout << "start: " << start_pt.transpose() << "\n:end: " <<
  // end_pt.transpose()
  //      << endl;
  ros::Time time_1, time_2;
  time_1 = ros::Time::now();
  double vis_time = 0.0;
  double transit_time1 = 0.0;
  double transit_time2 = 0.0;
  double get_neighbor_time = 0.0;
  ros::Time t1, t2;

  start_point = start_pt;
  end_point = end_pt;
  start_vel = start_v;
  start_acc = start_a;

  Vector3i start_idx = coord2gridIndex(start_pt);
  Vector3i end_idx = coord2gridIndex(end_pt);

  GridNodePtr startPtr = GridNodeMap[start_idx(0)][start_idx(1)][start_idx(2)];
  GridNodePtr endPtr = GridNodeMap[end_idx(0)][end_idx(1)][end_idx(2)];

  startPtr->state.head(3) = start_pt;
  startPtr->state.segment(3, 3) = start_v;
  // startPtr->state.tail(3) = Eigen::Vector3d::Zero(3);
  endPtr->state.head(3) = end_pt;
  endPtr->state.segment(3, 3) = end_v;

  startPtr->gScore = 0;
  startPtr->fScore = lambda_heu * getKinoDynamicHeu(startPtr, endPtr);
  // startPtr->fScore = 0.0;
  startPtr->id = 1;  // put start node in open set
  startPtr->coord = start_pt;

  openSet.clear();
  openSet.insert(
      make_pair(startPtr->fScore, startPtr));  // put start in open set

  GridNodePtr neighbor_ptr = NULL;
  GridNodePtr current_node = NULL;
  terminate_ptr = NULL;

  double tentative_gScore;
  int num_iter = 0;

  // num_ope = 0;
  // time_in_forward = 0.0;
  is_shot_succ = false;
  t_shot = 0.0;
  N_max_shot = 10;
  cnt_shot = 0;
  dis_shot = (startPtr->state.head(3) - endPtr->state.head(3)).norm();

  while (!openSet.empty()) {
    ++num_iter;
    current_node = openSet.begin()->second;

    int difference = 5;

    bool near_end =
        abs(current_node->index(0) - endPtr->index(0)) <= difference &&
        abs(current_node->index(1) - endPtr->index(1)) <= difference &&
        abs(current_node->index(2) - endPtr->index(2)) <= difference;

    bool reach_horizon =
        (current_node->state.head(3) - start_pt).norm() >= horizon;

    if (reach_horizon || near_end) {
      if (reach_horizon) cout << "[hybrid Astar]: Reach horizon" << endl;
      if (near_end) cout << "[hybrid Astar]: near end." << endl;

      if (near_end) {
        if (!shotHeu(current_node, endPtr))
          cout << "shot fail!" << endl;
        else
          cout << "shot successful!" << endl;
      }

      cout << "A* iteration: " << num_iter << endl;
      gridPath = retrievePath(current_node);
      terminate_ptr = current_node;
      has_path = true;

      // time_2 = ros::Time::now();
      // cout << "vis time:" << vis_time << endl;
      // cout << "Time consume in hybrid A star path finding is: "
      //      << (time_2 - time_1).toSec() - vis_time << endl;

      // cout << "time in outer state transition:" << transit_time1 << endl;
      // cout << "time in iner state transition:" << transit_time2 << endl;
      // cout << "time in get neighbor:" << get_neighbor_time << endl;
      // cout << "time in compare neighbor:" << compare_time << endl;

      // cout << "[hybrid Astar]: end-----------" << endl;
      if (near_end) {
        if (terminate_ptr->cameFrom == NULL && !is_shot_succ)
          return NO_PATH;
        else
          return REACH_END;
      } else if (reach_horizon)
        return REACH_HORIZON;
    }

    openSet.erase(openSet.begin());
    current_node->id = -1;  // move node from open to closed set.
    expandedNodes.push_back(current_node);

    // get the state of current node
    Eigen::VectorXd state = current_node->state;
    // cout << "\ncurernt state:" << state.transpose() << endl;

    // get neighbor of this node
    // t1 = ros::Time::now();

    Neighbors neighbors;
    if (init) {
      getNeighborInit(current_node, endPtr, neighbors, transit_time1,
                      transit_time2, vis_time);
      cout << "init num: " << neighbors.getData().size() << endl;
      init = false;
    } else {
      getNeighbor(current_node, endPtr, neighbors, transit_time1, transit_time2,
                  vis_time);
    }

    // t2 = ros::Time::now();
    // get_neighbor_time += (t2 - t1).toSec();

    vector<pair<Vector3i, KinoState>> neighbor_data = neighbors.getData();

    // iterate the neighbors
    for (int i = 0; i < int(neighbor_data.size()); ++i) {
      // get the neighbor node
      Eigen::Vector3i diff, neighbor_idx;
      diff = neighbor_data[i].first;
      neighbor_idx(0) = current_node->index(0) + diff(0);
      neighbor_idx(1) = current_node->index(1) + diff(1);
      neighbor_idx(2) = current_node->index(2) + diff(2);

      KinoState neighbor = neighbor_data[i].second;
      neighbor_ptr =
          GridNodeMap[neighbor_idx(0)][neighbor_idx(1)][neighbor_idx(2)];

      double edge_cost = neighbor.edge_cost;
      double heu = neighbor.heu;
      double optimal_time = neighbor.optimal_time;
      double duration = neighbor.duration;
      Eigen::Vector3d input = neighbor.input;
      Eigen::VectorXd x1 = neighbor.state;

      tentative_gScore = current_node->gScore + edge_cost;
      if (neighbor_ptr->index ==
          current_node->index)  // in the same grid, need compare
      {
        double current_fscore = current_node->fScore;
        double tentative_fscore = current_node->gScore + edge_cost + heu;

        if (tentative_fscore < current_fscore + tie_breaker) {
          // neighborPtr -> input = u;
          neighbor_ptr->state = x1;
          neighbor_ptr->gScore = tentative_gScore;
          neighbor_ptr->fScore = tentative_fscore;
          neighbor_ptr->optimal_time = optimal_time;
          neighbor_ptr->duration = duration;
          neighbor_ptr->input = input;

          if (neighbor_ptr->id == 1) {
            if (num_iter < 2) ROS_WARN(" re-insert, take place");

            openSet.erase(neighbor_ptr->nodeMapIt);
          }

          neighbor_ptr->id = 1;
          neighbor_ptr->nodeMapIt = openSet.insert(make_pair(
              neighbor_ptr->fScore,
              neighbor_ptr));  // put neighbor in open set and record it.
        }
      } else if (neighbor_ptr->id != -1)  // in close set
      {
        if (neighbor_ptr->id != 1)  // discover a new node
        {
          neighbor_ptr->id = 1;
          neighbor_ptr->state = x1;
          neighbor_ptr->optimal_time = optimal_time;
          neighbor_ptr->duration = duration;
          neighbor_ptr->input = input;
          neighbor_ptr->cameFrom = current_node;
          neighbor_ptr->gScore = tentative_gScore;
          neighbor_ptr->fScore = neighbor_ptr->gScore + heu;
          // neighbor_ptr->fScore = neighbor_ptr->gScore + 0.0;
          neighbor_ptr->nodeMapIt = openSet.insert(make_pair(
              neighbor_ptr->fScore,
              neighbor_ptr));  // put neighbor in open set and record it.
          continue;
        } else if (tentative_gScore <=
                   neighbor_ptr->gScore)  // already in open set, need compare
                                          // and update
        {
          neighbor_ptr->state = x1;
          neighbor_ptr->optimal_time = optimal_time;
          neighbor_ptr->duration = duration;
          neighbor_ptr->input = input;
          neighbor_ptr->cameFrom = current_node;
          neighbor_ptr->gScore = tentative_gScore;
          // neighbor_ptr->fScore = neighbor_ptr->gScore + 0.0;
          neighbor_ptr->fScore = tentative_gScore + heu;
          openSet.erase(neighbor_ptr->nodeMapIt);
          neighbor_ptr->nodeMapIt = openSet.insert(make_pair(
              neighbor_ptr->fScore,
              neighbor_ptr));  // put neighbor in open set and record it.
        }
      }
    }
  }

  time_2 = ros::Time::now();
  ROS_WARN("Can't find path. Time consume in hybrid A star path finding is %f",
           (time_2 - time_1).toSec());
  cout << "total number of iteration used in hybrid Astar: " << num_iter
       << endl;
  return NO_PATH;
}

vector<Vector3d> HybridAStarPathFinder::getPath() {
  vector<Vector3d> path;

  for (auto ptr : gridPath) path.push_back(ptr->coord);

  reverse(path.begin(), path.end());
  return path;
}

std::vector<GridNodePtr> HybridAStarPathFinder::getPathNodes() {
  return gridPath;
}

void HybridAStarPathFinder::resetPath() { gridPath.clear(); }

// The state is a vector of {x, y, z, vx, vy, vz}
void HybridAStarPathFinder::stateTransit1(Eigen::VectorXd init_state,
                                          Eigen::VectorXd& final_state,
                                          Eigen::Vector3d um, double tau) {
  // build the state transition matrix phi(tau)
  for (int i = 0; i < 3; ++i) phi(i, i + 3) = tau;

  // cout << "phi:\n" << phi << endl;

  // The build the integral terms
  Eigen::VectorXd integral(6);
  integral.head(3) = 0.5 * pow(tau, 2) * um;
  integral.tail(3) = tau * um;
  // cout << "integral:\n" << integral << endl;

  final_state = phi * init_state + integral;
  // cout << "The final state is:\n" << final_state.transpose() << endl << endl;
}

/* step is 1,2,3,4,5 */
void HybridAStarPathFinder::stateTransitPrecompute(Eigen::VectorXd pre_state,
                                                   Eigen::VectorXd& final_state,
                                                   Eigen::Vector3d um, int step,
                                                   double delta_t) {
  Eigen::VectorXd integral(6);
  double time = double(step) * delta_t;
  integral.head(3) = 0.5 * pow(time, 2) * um;
  integral.tail(3) = time * um;

  final_state = pre_state + integral;
}

void HybridAStarPathFinder::jerkInputStateTransit(Eigen::VectorXd init_state,
                                                  Eigen::VectorXd& final_state,
                                                  Eigen::Vector3d um,
                                                  double tau) {
  Eigen::Vector3d x0, v0, a0;
  x0 = init_state.head(3);
  v0 = init_state.segment(3, 3);
  a0 = init_state.tail(3);

  // add tau in transit matrix
  for (int i = 0; i < 6; ++i) phi(i, i + 3) = tau;
  for (int i = 0; i < 3; ++i) phi(i, i + 6) = 0.5 * tau * tau;
  // std::cout << "phi:\n" << phi << std::endl;

  // integral
  Eigen::VectorXd integral(9);
  integral.head(3) = (1 / 6.0) * pow(tau, 3) * um;
  integral.segment(3, 3) = 0.5 * pow(tau, 2) * um;
  integral.tail(3) = tau * um;
  // std::cout << "integral:\n" << integral << std::endl;

  // cout << "init:" << init_state << endl;
  final_state = phi * init_state + integral;
  // cout << "final: " << final_state.transpose() << endl;
}

void HybridAStarPathFinder::getNeighbor(
    GridNodePtr current_node, GridNodePtr end_node, Neighbors& neighbors,
    double& transit_time1, double& transit_time2, double& vis_time) {
  ros::Time t1, t2;
  // applying different um and tau, we get a series of final state from current
  // state
  GridNodePtr nptr;
  Eigen::VectorXd state = current_node->state;

  Eigen::Vector3d um;
  double res = 1 / 2.0, time_res = 1 / 1.0;
  int pri_num = 0;

  /* precompute the state transition matrix at different time stamp */
  const int check_num = 10;
  vector<Eigen::VectorXd> precompute_states;
  precompute_states.resize(check_num);
  double delta_tau = (1 / double(check_num)) * max_tau;
  for (int j = 0; j < check_num; ++j) {
    double dt = double(j + 1) * delta_tau;
    Eigen::VectorXd pre_state(6);
    Eigen::MatrixXd pre_phi = Eigen::MatrixXd::Zero(6, 6);
    pre_phi.block(3, 0, 3, 3) = Eigen::MatrixXd::Zero(3, 3);
    pre_phi.block(0, 0, 3, 3) = pre_phi.block(3, 3, 3, 3) =
        Eigen::MatrixXd::Identity(3, 3);
    for (int i = 0; i < 3; ++i) pre_phi(i, i + 3) = dt;
    pre_state = pre_phi * state;
    precompute_states[j] = pre_state;
  }

  for (double ax = -max_acc; ax <= max_acc + 1e-3; ax += max_acc * res)
    for (double ay = -max_acc; ay <= max_acc + 1e-3; ay += max_acc * res)
      for (double az = -max_acc; az <= max_acc + 1e-3; az += max_acc * res) {
        um << ax, ay, 0.5 * az;
        // state transit
        Eigen::VectorXd x1;

        double out_transit_t = 0.0, if_t = 0.0, collision_t = 0.0, cal_t = 0.0,
               compare_t = 0.0, heu_t = 0.0;
        // t1 = ros::Time::now();
        // stateTransit1(state, x1, um, max_tau);
        stateTransitPrecompute(precompute_states[check_num - 1], x1, um,
                               check_num, delta_tau);

        // t2 = ros::Time::now();
        // transit_time1 += (t2 - t1).toSec();
        // out_transit_t = (t2 - t1).toSec();
        // cout << "outer state transit: " << (t2 - t1).toSec() * 1000 << endl;

        // jerkInputStateTransit(state, x1, um, max_tau);

        // /* // display the transit */
        // ros::Time t1 = ros::Time::now();
        // static int id_temp = 0;
        // vector<Eigen::Vector3d> path_temp;
        // for (double dt = 0; dt < max_tau + 1e-3; dt += 0.005) {
        //   Eigen::VectorXd xt;
        //   // jerkInputStateTransit(state, xt, um, dt);
        //   stateTransit1(state, xt, um, dt);
        //   path_temp.push_back(xt.head(3));
        // }
        // displayPathWithColor(path_temp, 0.01, 1, id_temp % 2000);
        // ++id_temp;
        // ros::Time t2 = ros::Time::now();
        // vis_time += (t2 - t1).toSec();

        // cout << "state:" << x1.transpose() << endl;

        // stay in local range

        // t1 = ros::Time::now();

        Eigen::Vector3i idx1 = coord2gridIndex(x1.head(3));
        if (idx1(0) < 0 || idx1(0) >= grid_num(0) || idx1(1) < 0 ||
            idx1(1) >= grid_num(1) || idx1(2) < 0 || idx1(2) >= grid_num(2)) {
          // cout << "not in range" << endl;
          break;
        }

        // vel feasible
        Eigen::Vector3d v1 = x1.segment(3, 3);
        if (fabs(v1(0)) > max_vel || fabs(v1(1)) > max_vel ||
            fabs(v1(2)) > max_vel) {
          // cout << "vel end not feasible" << endl;
          break;
        }

        // check if it is neighbor
        Eigen::Vector3i diff = idx1 - current_node->index;
        if (diff.norm() == 0) continue;
        // cout << "neighbor:" << diff.transpose() << endl;

        // t2 = ros::Time::now();
        // cout << "if time: " << (t2 - t1).toSec() * 1000 << endl;
        // if_t = (t2 - t1).toSec();

        // ros::Time t3 = ros::Time::now();
        // collision free
        Eigen::Vector3i idt;
        Eigen::VectorXd xt;

        idt = coord2gridIndex(x1.head(3));
        nptr = GridNodeMap[idt(0)][idt(1)][idt(2)];
        if (nptr->occupancy > 0.5) {
          // cout << "collision" << endl;
          break;
        }

        bool is_occ = false;
        for (int j = 0; j < check_num - 1; ++j) {
          double dt = double(j + 1) * delta_tau;

          // t1 = ros::Time::now();
          // stateTransit1(state, xt, um, dt);
          stateTransitPrecompute(precompute_states[j], xt, um, j + 1,
                                 delta_tau);
          // t2 = ros::Time::now();
          // collision_t = (t2 - t1).toSec();
          // transit_time2 += (t2 - t1).toSec();

          idt = coord2gridIndex(xt.head(3));
          nptr = GridNodeMap[idt(0)][idt(1)][idt(2)];
          // nptr->occupancy > 0.5
          if (nptr->distance <= 0.2) {
            is_occ = true;
            // cout << "collision, " << nptr->distance << endl;
            break;
          }
        }
        if (is_occ) break;
        // ros::Time t4 = ros::Time::now();
        // cout << "check collision: " << (t4 - t3).toSec() * 1000 << endl;
        // collision_t = (t4 - t3).toSec();

        // caluculate f_score
        // t1 = ros::Time::now();

        double optimal_time;
        KinoState candidate;
        candidate.edge_cost = (um.squaredNorm() + w_time) * max_tau;

        // t3 = ros::Time::now();
        candidate.heu =
            lambda_heu * getKinoDynamicHeu(x1, end_node->state, optimal_time);
        // t4 = ros::Time::now();
        // heu_t = (t4 - t3).toSec();

        candidate.optimal_time = optimal_time;
        candidate.state = x1;
        candidate.input = um;
        candidate.duration = max_tau;

        // cout << "\ninput:" << um.transpose() << ", edge cost:" <<
        // candidate.edge_cost
        //      << ", heu:" << candidate.heu << ", sum cost:" <<
        //      candidate.edge_cost + candidate.heu << "\n";
        // t2 = ros::Time::now();
        // cout << "cal something: " << (t2 - t1).toSec() * 1000 << endl;
        // cal_t = (t2 - t1).toSec();

        // t1 = ros::Time::now();
        KinoState exist_neighbor;
        // if not found, insert it directly
        if (!neighbors.find(diff, exist_neighbor)) {
          neighbors.add(diff, candidate);
          // cout << "add new " << endl;
        } else  // compare the edge_cost + heu_cost
        {
          bool replace = (candidate.edge_cost + candidate.heu) <
                         (exist_neighbor.edge_cost + exist_neighbor.heu);
          if (replace) {
            neighbors.erase(diff);
            neighbors.add(diff, candidate);
            // cout << "replace" << endl;
          }
        }

        // t2 = ros::Time::now();
        // compare_time += (t2 - t1).toSec();
        // compare_t = (t2 - t1).toSec();

        /* statistic time */
        // double tot_t = out_transit_t + if_t + collision_t + cal_t +
        // compare_t; out_transit_t /= tot_t; if_t /= tot_t; collision_t /=
        // tot_t; cal_t /= tot_t; compare_t /= tot_t; heu_t /= tot_t; double sum
        // = out_transit_t + if_t + collision_t + cal_t + compare_t;

        // cout << fixed << setprecision(3) << "\n\nout transit: " <<
        // out_transit_t
        //      << "\nif_t: " << if_t << "\ncollision: " << collision_t
        //      << "\ncal: " << cal_t << "\nheu t: " << heu_t
        //      << "\ncompare: " << compare_t << "\nsum: " << sum << endl;
      }
}

double HybridAStarPathFinder::getKinoDynamicHeu(GridNodePtr node1,
                                                GridNodePtr node2) {
  ros::Time t1, t2;

  // t1 = ros::Time::now();

  const Vector3d dp = node2->state.head(3) - node1->state.head(3);
  // const Vector3d v0 = node1->state.segment(3, 3);
  const Vector3d v0 = node1->state.segment(3, 3);
  const Vector3d v1 = node2->state.segment(3, 3);

  double c1 = -36 * dp.dot(dp);
  double c2 = 24 * (v0 + v1).dot(dp);
  double c3 = -4 * (v0.dot(v0) + v0.dot(v1) + v1.dot(v1));
  double c4 = 0;
  double c5 = w_time;
  // t2 = ros::Time::now();
  // cout << "1: " << (t2 - t1).toSec() * 1000 << endl;

  // t1 = ros::Time::now();

  std::vector<double> ts = quartic(c5, c4, c3, c2, c1);
  /* try not consider w_time */
  // std::vector<double> ts = quadratic(c1, c2, c3);
  // t2 = ros::Time::now();
  // cout << "2: " << (t2 - t1).toSec() * 1000 << endl;

  // t1 = ros::Time::now();

  double v_max = max_vel;
  double t_bar =
      (node1->state.head(3) - node2->state.head(3)).lpNorm<Infinity>() / v_max;

  ts.push_back(t_bar);
  // t2 = ros::Time::now();
  // cout << "3: " << (t2 - t1).toSec() * 1000 << endl;

  // t1 = ros::Time::now();

  double cost = std::numeric_limits<double>::max();
  double t_d = t_bar;

  for (auto t : ts) {
    if (t < t_bar) continue;
    double c = -c1 / (3 * t * t * t) - c2 / (2 * t * t) - c3 / t + w_time * t;
    if (c < cost) {
      cost = c;
      t_d = t;
    }
  }

  // node1->optimal_t = t_d;
  node1->optimal_time = t_d;

  // t2 = ros::Time::now();
  // cout << "4: " << (t2 - t1).toSec() * 1000 << endl;

  return 1.0 * (1 + tie_breaker) * cost;
}

double HybridAStarPathFinder::getKinoDynamicHeu(Eigen::VectorXd x1,
                                                Eigen::VectorXd x2,
                                                double& optimal_time) {
  const Vector3d dp = x2.head(3) - x1.head(3);
  const Vector3d v0 = x1.segment(3, 3);
  const Vector3d v1 = x2.segment(3, 3);

  double c1 = -36 * dp.dot(dp);
  double c2 = 24 * (v0 + v1).dot(dp);
  double c3 = -4 * (v0.dot(v0) + v0.dot(v1) + v1.dot(v1));
  double c4 = 0;
  double c5 = w_time;

  std::vector<double> ts = quartic(c5, c4, c3, c2, c1);

  double v_max = max_vel;
  // double t_bar = (x1.head(3) - x2.head(3)).lpNorm<Infinity>() / v_max;
  Eigen::Vector3d to_goal = dp / v_max;
  double t_bar = max(fabs(to_goal(0)), fabs(to_goal(1)));
  t_bar = max(t_bar, fabs(to_goal(2)));
  ts.push_back(t_bar);

  double cost = 100000000;
  double t_d = t_bar;

  for (auto t : ts) {
    if (t < t_bar) continue;
    double c = -c1 / (3 * t * t * t) - c2 / (2 * t * t) - c3 / t + w_time * t;
    if (c < cost) {
      cost = c;
      t_d = t;
    }
  }

  // node1->optimal_t = t_d;
  optimal_time = t_d;

  return 1.0 * (1 + tie_breaker) * cost;
}

vector<double> HybridAStarPathFinder::cubic(double a, double b, double c,
                                            double d) {
  vector<double> dts;

  double a2 = b / a;
  double a1 = c / a;
  double a0 = d / a;

  double Q = (3 * a1 - a2 * a2) / 9;
  double R = (9 * a1 * a2 - 27 * a0 - 2 * a2 * a2 * a2) / 54;
  double D = Q * Q * Q + R * R;
  if (D > 0) {
    double S = std::cbrt(R + sqrt(D));
    double T = std::cbrt(R - sqrt(D));
    dts.push_back(-a2 / 3 + (S + T));
    return dts;
  } else if (D == 0) {
    double S = std::cbrt(R);
    dts.push_back(-a2 / 3 + S + S);
    dts.push_back(-a2 / 3 - S);
    return dts;
  } else {
    double theta = acos(R / sqrt(-Q * Q * Q));
    dts.push_back(2 * sqrt(-Q) * cos(theta / 3) - a2 / 3);
    dts.push_back(2 * sqrt(-Q) * cos((theta + 2 * M_PI) / 3) - a2 / 3);
    dts.push_back(2 * sqrt(-Q) * cos((theta + 4 * M_PI) / 3) - a2 / 3);
    return dts;
  }
}

vector<double> HybridAStarPathFinder::quartic(double a, double b, double c,
                                              double d, double e) {
  vector<double> dts;

  double a3 = b / a;
  double a2 = c / a;
  double a1 = d / a;
  double a0 = e / a;

  vector<double> ys =
      cubic(1, -a2, a1 * a3 - 4 * a0, 4 * a2 * a0 - a1 * a1 - a3 * a3 * a0);
  double y1 = ys.front();
  double r = a3 * a3 / 4 - a2 + y1;
  if (r < 0) return dts;

  double R = sqrt(r);
  double D, E;
  if (R != 0) {
    D = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 +
             0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
    E = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 -
             0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
  } else {
    D = sqrt(0.75 * a3 * a3 - 2 * a2 + 2 * sqrt(y1 * y1 - 4 * a0));
    E = sqrt(0.75 * a3 * a3 - 2 * a2 - 2 * sqrt(y1 * y1 - 4 * a0));
  }

  if (!std::isnan(D)) {
    dts.push_back(-a3 / 4 + R / 2 + D / 2);
    dts.push_back(-a3 / 4 + R / 2 - D / 2);
  }
  if (!std::isnan(E)) {
    dts.push_back(-a3 / 4 - R / 2 + E / 2);
    dts.push_back(-a3 / 4 - R / 2 - E / 2);
  }

  return dts;
}

vector<double> HybridAStarPathFinder::quadratic(double a, double b, double c) {
  vector<double> roots;

  /* delta */
  double delta = b * b - 4 * a * c;

  if (delta < 0)
    ;
  else if (delta == 0) {
    roots.push_back(-b / (2 * a));
  } else if (delta > 0) {
    roots.push_back((-b + sqrt(delta)) / (2 * a));
    roots.push_back((-b - sqrt(delta)) / (2 * a));
  }
  return roots;
}

double HybridAStarPathFinder::getOptimalTime(Eigen::Vector3d p0,
                                             Eigen::Vector3d p1,
                                             Eigen::Vector3d v0) {
  Eigen::Vector3d dp = p1 - p0;

  /* solve candidate optimal time */
  double a = 3 * v0.dot(v0);
  double b = -12 * dp.dot(v0);
  double c = 9 * dp.dot(dp);
  vector<double> roots = quadratic(a, b, c);

  /* select min cost time */
  double min_cost = 1000000;
  double opti_time;
  for (int i = 0; i < int(roots.size()); ++i) {
    if (roots[i] < 0) continue;

    double T = roots[i];
    double cost = 3 * (v0 * T - dp).dot(v0 * T - dp) / (T * T * T);
    if (cost < min_cost) {
      min_cost = cost;
      opti_time = T;
    }
  }

  return opti_time;
}

/* get min acc integral one shot; check feasibility and reset time */
Eigen::MatrixXd HybridAStarPathFinder::getShotTrajectory(Eigen::Vector3d p0,
                                                         Eigen::Vector3d p1,
                                                         Eigen::Vector3d v0,
                                                         Eigen::Vector3d& v1,
                                                         double& T) {
  Eigen::Vector3d dp = p1 - p0;

  /* cal initial optimal time */
  // double T = 1.5 * dp.norm() / max_vel;
  T = getOptimalTime(p0, p1, v0);

  /* check feasibility and reset time */
  Eigen::Vector3d ve = v0 + 3 * (dp - v0 * T) / (2 * T);
  for (int i = 0; i < 3; ++i) {
    if (ve(i) > (2.5 / 3) * max_vel) {
      double Tp = 3 * dp(i) / (2 * (max_vel + 0.5 * v0(i)));
      if (Tp > T) T = Tp;
    }
  }
  v1 = v0 + 3 * (dp - v0 * T) / (2 * T);

  /* get coefficient */
  Eigen::Vector3d a = -(dp - v0 * T) / 2.0 / T / T / T;
  Eigen::Vector3d b = 3 * (dp - v0 * T) / 2.0 / T / T;
  Eigen::Vector3d c = v0;
  Eigen::Vector3d d = p0;
  // SETCY << "a:" << a.transpose() << ", b:" << b.transpose()
  //       << ", c:" << c.transpose() << ", d" << d.transpose() << endl;

  Eigen::MatrixXd coef(3, 4);
  coef.col(3) = a;
  coef.col(2) = b;
  coef.col(1) = c;
  coef.col(0) = d;

  return coef;
}

bool HybridAStarPathFinder::freeEndVelShot(GridNodePtr node1,
                                           GridNodePtr node2) {
  const Vector3d p0 = node1->state.head(3);
  const Vector3d p1 = node2->state.head(3);
  const Vector3d v0 = node1->state.segment(3, 3);

  /* get shot trajectory and its optimal time, this shot is always kinodynamic
   * feasible!! */
  double T;
  Eigen::Vector3d v1;
  Eigen::MatrixXd coef = getShotTrajectory(p0, p1, v0, v1, T);
  end_vel = v1;

  /* check feasibility */
  Vector3d coord;
  VectorXd poly1d, t;
  Vector3i index;

  /* get a proper delta_t */
  int step = ceil((p1 - p0).norm() / (sqrt(3) * resolution));
  double t_delta = T / double(step);

  for (double time = t_delta; time <= T; time += t_delta) {
    t = VectorXd::Zero(4);
    for (int j = 0; j < 4; j++) t(j) = pow(time, j);

    for (int dim = 0; dim < 3; dim++) {
      poly1d = coef.row(dim);
      coord(dim) = poly1d.dot(t);
    }
    index = coord2gridIndex(coord);

    if (index(0) < 0 || index(0) >= grid_num(0) || index(1) < 0 ||
        index(1) >= grid_num(1) || index(2) < 0 || index(2) >= grid_num(2)) {
      cout << "shot outside map" << endl;
      return false;
    }

    if (GridNodeMap[index(0)][index(1)][index(2)]->occupancy >
        0.5)  // collision
    {
      cout << "shot in obstacles" << endl;
      return false;
    }
  }

  coef_shot = coef;
  t_shot = T;
  is_shot_succ = true;

  cout << "shot success, end vel: " << end_vel.transpose() << endl;
  return true;
}

bool HybridAStarPathFinder::shotHeu(GridNodePtr node1, GridNodePtr node2) {
  // cnt_shot++;
  // if (cnt_shot < N_max_shot)
  //     return false;
  // else
  // {
  //     cnt_shot = 0;
  //     N_max_shot = ceil((node1->state.head(3) - node2->state.head(3)).norm()
  //     / dis_shot * N_max);
  // }
  // TODO: need to check the velocity and acceleration feasibility
  //  ****** now check the feasibility of the optimal polynomial, by using t_d

  const Vector3d p0 = node1->state.head(3);
  const Vector3d dp = node2->state.head(3) - p0;
  const Vector3d v0 = node1->state.segment(3, 3);
  const Vector3d v1 = node2->state.segment(3, 3);
  const Vector3d dv = v1 - v0;
  double t_d = node1->optimal_time;
  MatrixXd coef(3, 4);
  end_vel = v1;

  /* get coefficient */
  Vector3d a =
      1.0 / 6.0 *
      (-12.0 / (t_d * t_d * t_d) * (dp - v0 * t_d) + 6 / (t_d * t_d) * dv);
  Vector3d b = 0.5 * (6.0 / (t_d * t_d) * (dp - v0 * t_d) - 2 / t_d * dv);
  Vector3d c = v0;
  Vector3d d = p0;

  coef.col(3) = a;
  coef.col(2) = b;
  coef.col(1) = c;
  coef.col(0) = d;

  // *** the OPTIMAL polynomial is : 1/6 * alpha * t^3 + 1/2 * beta * t^2 + v0
  // * t + p0; denote as : a*t^3 + b*t^2
  // + v0*t + p0
  Vector3d coord, vel, acc;
  VectorXd poly1d, t, polyv, polya;
  Vector3i index;

  Eigen::MatrixXd Tm(4, 4);
  Tm << 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0;

  // forward checking of trajectory
  double t_delta = t_d / 10;
  for (double time = t_delta; time <= t_d; time += t_delta) {
    t = VectorXd::Zero(4);
    for (int j = 0; j < 4; j++) t(j) = pow(time, j);

    for (int dim = 0; dim < 3; dim++) {
      poly1d = coef.row(dim);
      coord(dim) = poly1d.dot(t);
      vel(dim) = (Tm * poly1d).dot(t);
      acc(dim) = (Tm * Tm * poly1d).dot(t);

      if (fabs(vel(dim)) > max_vel || fabs(acc(dim)) > max_acc) {
        // cout << "vel:" << vel(dim) << ", acc:" << acc(dim) << endl;
      }
    }

    index = coord2gridIndex(coord);

    if (index(0) < 0 || index(0) >= grid_num(0) || index(1) < 0 ||
        index(1) >= grid_num(1) || index(2) < 0 || index(2) >= grid_num(2)) {
      return false;
    }

    if (GridNodeMap[index(0)][index(1)][index(2)]->occupancy > 0.5) {
      return false;
    }
  }
  coef_shot = coef;
  t_shot = t_d;
  is_shot_succ = true;
  // cout << "shot success!" << endl;
  return true;
}

void HybridAStarPathFinder::getKinoTrajMat(Eigen::MatrixXd& Pos,
                                           Eigen::MatrixXd& Vel,
                                           Eigen::MatrixXd& Acc,
                                           Eigen::VectorXd& Time) {
  vector<Eigen::Vector3d> vec_pos, vec_vel, vec_acc;
  vector<double> vec_time;

  GridNodePtr ptr = terminate_ptr;
  if (ptr == NULL) return;  // ok not to check one shot in benckmark

  /* add shot first */
  if (is_shot_succ) {
    vec_time.push_back(t_shot);
    vec_pos.push_back(end_point);
    vec_vel.push_back(Eigen::Vector3d(0, 0, 0));
    vec_acc.push_back(Eigen::Vector3d(0, 0, 0));

    vec_pos.push_back(ptr->state.head(3));
    vec_vel.push_back(ptr->state.tail(3));
    vec_acc.push_back(Eigen::Vector3d(0, 0, 0));
  } else {
    vec_time.push_back(ptr->duration);
    vec_pos.push_back(ptr->state.head(3));
    vec_vel.push_back(ptr->state.tail(3));
    vec_acc.push_back(ptr->input);

    vec_pos.push_back(ptr->cameFrom->state.head(3));
    vec_vel.push_back(ptr->cameFrom->state.tail(3));
    vec_acc.push_back(ptr->input);
    ptr = ptr->cameFrom;
  }

  while (ptr->cameFrom != NULL) {
    vec_time.push_back(ptr->duration);
    vec_pos.push_back(ptr->cameFrom->state.head(3));
    vec_vel.push_back(ptr->cameFrom->state.tail(3));
    vec_acc.push_back(ptr->input);
    ptr = ptr->cameFrom;
  }

  reverse(vec_time.begin(), vec_time.end());
  reverse(vec_pos.begin(), vec_pos.end());
  reverse(vec_vel.begin(), vec_vel.end());
  reverse(vec_acc.begin(), vec_acc.end());

  /* assign to matrix */
  Time.resize(vec_time.size());
  Pos.resize(vec_pos.size(), 3);
  Vel.resize(vec_vel.size(), 3);
  Acc.resize(vec_acc.size(), 3);

  for (int i = 0; i < vec_time.size(); ++i) Time(i) = vec_time[i];

  for (int i = 0; i < vec_pos.size(); ++i) {
    Pos.row(i) = vec_pos[i];
    Vel.row(i) = vec_vel[i];
    Acc.row(i) = vec_acc[i];
  }
}

vector<Eigen::Vector3d> HybridAStarPathFinder::getKinoTraj(double delta_t) {
  vector<Vector3d> state_list;

  GridNodePtr ptr = terminate_ptr;
  if (ptr == NULL)  // no path found
    return state_list;

  VectorXd xt, xu;

  // ROS_WARN("[hybridAstarSearch] check point's index");
  while (ptr->cameFrom != NULL) {
    Vector3d u = ptr->input;
    double duration = ptr->duration;
    // cout << "index:" << ptr->index.transpose() << ", duration:" << duration
    // << endl;

    xt = ptr->cameFrom->state;
    for (double t = duration; t >= -1e-3; t -= delta_t) {
      stateTransit1(xt, xu, u, t);
      // jerkInputStateTransit(xt, xu,gg u, t);
      // cout << "t:" << t << ", state:" << xu.head(3).transpose() << endl;
      state_list.push_back(xu.head(3));
    }

    ptr = ptr->cameFrom;
  }

  reverse(state_list.begin(), state_list.end());

  if (is_shot_succ)  // add shoting heuristic trajectory to the kino traj list
  {
    Vector3d coord;
    VectorXd poly1d, t;

    for (double time = delta_t; time <= t_shot; time += delta_t) {
      for (int dim = 0; dim < 3; dim++) {
        poly1d = coef_shot.row(dim);
        t = VectorXd::Zero(4);

        for (int j = 0; j < 4; j++) t(j) = pow(time, j);

        coord(dim) = poly1d.dot(t);
      }

      state_list.push_back(coord);
    }
  }

  return state_list;
}

// This function return samples in uniform time interval
// input: ts', N(sample num)
// output: ts,  mat of 3 x (K+1)*(N+1) and each row for t, x, y, z, 3 x 3
// matrix start_end_acc
Eigen::MatrixXd HybridAStarPathFinder::getSamples(double& ts, int& K, int N,
                                                  bool repeat) {
  Eigen::MatrixXd samples;

  GridNodePtr ptr = terminate_ptr;
  if (ptr == NULL) {
    cout << "no path found, return null sample" << endl;
    return samples;
  }

  // cal the accumulated time of the path
  double T = 0.0;
  if (is_shot_succ) T += t_shot;

  // cout << "is shot success: " << is_shot_succ << "， t shot: " << t_shot
  //      << endl;

  while (ptr->cameFrom != NULL) {
    T += ptr->duration;
    ptr = ptr->cameFrom;
  }

  // cout << "accumulated time:" << T << endl;

  // cal ts, tm
  K = floor(T / ts);
  ts = T / (K + 1);
  double tm = ts / N;

  // cout << "K:" << K << ", N:" << N << ", ts:" << ts << ", tm:" << tm << endl;

  // get samples
  // bool in_shot = true;
  bool in_shot = is_shot_succ;
  int Nt = 0;
  Eigen::VectorXd sx, sy, sz;
  if (repeat) {
    sx.resize((N + 1) * (K + 1));
    sy.resize((N + 1) * (K + 1));
    sz.resize((N + 1) * (K + 1));
  } else {
    sx.resize(K + 2);
    sy.resize(K + 2);
    sz.resize(K + 2);
  }
  double T_accumulate = T, t;
  if (is_shot_succ)
    t = t_shot;
  else {
    t = terminate_ptr->duration;
    ptr = terminate_ptr;
    end_vel = ptr->state.tail(3);
  }

  while (true) {
    if (in_shot)  // cal one shot coordinate
    {
      Vector3d coord;
      VectorXd poly1d, tv;
      for (int dim = 0; dim < 3; dim++) {
        poly1d = coef_shot.row(dim);
        tv = VectorXd::Zero(4);

        for (int j = 0; j < 4; j++) tv(j) = pow(t, j);

        coord(dim) = poly1d.dot(tv);
      }

      sx(Nt) = coord(0);
      sy(Nt) = coord(1);
      sz(Nt) = coord(2);
      ++Nt;
      // segment connecting point must be added twice
      if (repeat && Nt % (N + 1) == 0 && Nt != (K + 1) * (N + 1)) {
        sx(Nt) = coord(0);
        sy(Nt) = coord(1);
        sz(Nt) = coord(2);
        ++Nt;
      }

      // move to next sample
      t -= tm;
      T_accumulate -= tm;
      if (t < -1e-5)  // outside the range of oneshot
      {
        in_shot = false;
        ptr = terminate_ptr;

        /* only shot, no search */
        if (ptr->cameFrom == NULL) {
          samples.resize(3, K + 5);
          samples.block(0, 0, 1, K + 2) = sx.reverse().transpose();
          samples.block(1, 0, 1, K + 2) = sy.reverse().transpose();
          samples.block(2, 0, 1, K + 2) = sz.reverse().transpose();
          samples.col(K + 2) = start_vel;
          samples.col(K + 3) = end_vel;
          samples.col(K + 4) = 2.0 * coef_shot.col(2);
          // cout << "return from shot" << endl;
          return samples;
        }

        t += ptr->duration;
        // cout << "go to search part" << endl;
      }
    } else  // cal coordinate of normal path
    {
      Eigen::VectorXd xu;
      Vector3d u = ptr->input;
      Eigen::VectorXd xt = ptr->cameFrom->state;
      stateTransit1(xt, xu, u, t);

      sx(Nt) = xu(0);
      sy(Nt) = xu(1);
      sz(Nt) = xu(2);
      ++Nt;

      // segment connecting point must be added twice
      if (repeat && Nt % (N + 1) == 0 && Nt != (K + 1) * (N + 1)) {
        sx(Nt) = xu(0);
        sy(Nt) = xu(1);
        sz(Nt) = xu(2);
        ++Nt;
      }

      // move to next sample
      t -= tm;
      T_accumulate -= tm;

      if (t < -1e-5)  // outside the range of path segment
      {
        if (ptr->cameFrom->cameFrom ==
            NULL)  // reach the first node, finish all samples
        {
          if (repeat)
            samples.resize(3, (K + 1) * (N + 1));
          else
            samples.resize(3, K + 5);
          samples.block(0, 0, 1, K + 2) = sx.reverse().transpose();
          samples.block(1, 0, 1, K + 2) = sy.reverse().transpose();
          samples.block(2, 0, 1, K + 2) = sz.reverse().transpose();
          samples.col(K + 2) = start_vel;
          samples.col(K + 3) = end_vel;
          samples.col(K + 4) = ptr->input;

          // get start acc
          return samples;
        } else  // not reach the first node
        {
          ptr = ptr->cameFrom;
          t += ptr->duration;
          // cout << "input:" << ptr->input.transpose() << ", velocity" <<
          // ptr->state.tail(3).transpose()
          //      << endl;
        }
      }
    }
  }
}

void HybridAStarPathFinder::getNeighborInit(
    GridNodePtr current_node, GridNodePtr end_node, Neighbors& neighbors,
    double& transit_time1, double& transit_time2, double& vis_time) {
  ros::Time t1, t2;
  GridNodePtr nptr;
  Eigen::VectorXd state = current_node->state;

  Eigen::Vector3d um;
  double res = 1 / 8.0;
  int pri_num = 0;

  for (double tau = res * init_max_tau; tau <= init_max_tau;
       tau += res * init_max_tau) {
    um = start_acc;
    Eigen::VectorXd x1;
    stateTransit1(state, x1, um, tau);

    Eigen::Vector3i idx1 = coord2gridIndex(x1.head(3));
    if (idx1(0) < 0 || idx1(0) >= grid_num(0) || idx1(1) < 0 ||
        idx1(1) >= grid_num(1) || idx1(2) < 0 || idx1(2) >= grid_num(2)) {
      continue;
    }

    // vel feasible
    Eigen::Vector3d v1 = x1.segment(3, 3);
    if (fabs(v1(0)) > max_vel || fabs(v1(1)) > max_vel || fabs(v1(2)) > max_vel)
      continue;

    // check if it is neighbor
    Eigen::Vector3i diff = idx1 - current_node->index;
    if (diff.norm() == 0) continue;

    // collision free
    Eigen::Vector3i idt;
    Eigen::VectorXd xt;
    bool is_occ = false;
    for (double tck = tau / 10.0; tck <= tau; tck += tau / 10.0) {
      stateTransit1(state, xt, um, tck);
      idt = coord2gridIndex(xt.head(3));
      nptr = GridNodeMap[idt(0)][idt(1)][idt(2)];
      if (nptr->distance <= 0.2) {
        is_occ = true;
        break;
      }
    }
    if (is_occ) continue;

    /* add to neighbor vector */
    double optimal_time;
    KinoState candidate;
    candidate.edge_cost = (um.squaredNorm() + w_time) * tau;
    candidate.heu =
        lambda_heu * getKinoDynamicHeu(x1, end_node->state, optimal_time);
    candidate.optimal_time = optimal_time;
    candidate.state = x1;
    candidate.input = um;
    candidate.duration = tau;

    KinoState exist_neighbor;
    if (!neighbors.find(diff, exist_neighbor)) {
      neighbors.add(diff, candidate);
    } else {
      bool replace = (candidate.edge_cost + candidate.heu) <
                     (exist_neighbor.edge_cost + exist_neighbor.heu);
      if (replace) {
        neighbors.erase(diff);
        neighbors.add(diff, candidate);
      }
    }
  }
}