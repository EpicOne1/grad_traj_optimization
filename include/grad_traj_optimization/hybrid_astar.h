#ifndef _HYBRID_ASTAR_H
#define _HYBRID_ASTAR_H

#include <ros/console.h>
#include <ros/ros.h>
#include <Eigen/Eigen>
#include <iostream>
#include <map>
#include <string>
#include "grad_traj_optimization/data_type.h"
#include "grad_traj_optimization/sdf_map.h"

#define REACH_HORIZON 1
#define REACH_END 2
#define NO_PATH 3

class HybridAStarPathFinder {
 private:
  struct KinoState {
    double edge_cost;
    double heu;
    double optimal_time;
    double duration;
    Eigen::VectorXd input;
    Eigen::VectorXd state;
  };

  class Neighbors {
   private:
    /* data */
    std::vector<std::pair<Eigen::Vector3i, KinoState>> neighbor_data;

   public:
    Neighbors(/* args */) {}
    ~Neighbors() {}

    void erase(Eigen::Vector3i id) {
      for (int i = 0; i < int(neighbor_data.size()); ++i) {
        if ((id - neighbor_data[i].first).norm() == 0) {
          neighbor_data.erase(neighbor_data.begin() + i);
          break;
        }
      }
    }

    void add(Eigen::Vector3i id, KinoState state) {
      neighbor_data.push_back(std::make_pair(id, state));
    }

    bool find(Eigen::Vector3i id, KinoState& state) {
      for (int i = 0; i < int(neighbor_data.size()); ++i) {
        if ((id - neighbor_data[i].first).norm() == 0) {
          state = neighbor_data[i].second;
          return true;
        }
      }

      return false;
    }

    std::vector<std::pair<Eigen::Vector3i, KinoState>> getData() {
      return this->neighbor_data;
    }
  };

  /* statistic */
  double compare_time = 0.0;

  /* kinodynamic */
  Eigen::MatrixXd phi;    // state transit matrix
  double max_tau = 0.25;  // transition time
  double init_max_tau = 0.8;
  double max_vel = 3.0;
  double max_acc = 3.0;
  double max_jerk = 3.0;
  double w_time = 10.0;
  Eigen::Vector3d start_point, end_point, start_vel, end_vel, start_acc,
      end_acc;
  double horizon;
  double lambda_heu;

  enum INPUT_TYPE { JERK_INPUT, ACC_INPUT };
  int input_type;
  int visualization;

  /* graph search */
  std::vector<GridNodePtr> expandedNodes;
  std::vector<GridNodePtr> gridPath;
  GridNodePtr*** GridNodeMap;
  std::multimap<double, GridNodePtr> openSet;

  /* map */
  double resolution, inv_resolution;
  double tie_breaker = 1.0 + 1.0 / 10000;
  Eigen::Vector3i grid_num;
  Eigen::Vector3d origin;

  Eigen::Vector3d gridIndex2coord(Eigen::Vector3i index);
  Eigen::Vector3i coord2gridIndex(Eigen::Vector3d pt);
  GridNodePtr pos2gridNodePtr(Eigen::Vector3d pos);
  GridNodePtr terminate_ptr;

  /* one shot */
  bool is_shot_succ = false;
  int N_max_shot = 10;
  int N_max = 10;
  int cnt_shot = 0;
  Eigen::MatrixXd coef_shot;
  double t_shot, dis_shot;
  bool has_path = false;

  /* old heuristic */
  double getDiagHeu(GridNodePtr node1, GridNodePtr node2);
  double getManhHeu(GridNodePtr node1, GridNodePtr node2);
  double getEuclHeu(GridNodePtr node1, GridNodePtr node2);
  double getHeu(GridNodePtr node1, GridNodePtr node2);
  bool shotHeu(GridNodePtr node1, GridNodePtr node2);
  double getKinoDynamicHeu(GridNodePtr node1, GridNodePtr node2);
  double getKinoDynamicHeu(Eigen::VectorXd x1, Eigen::VectorXd x2,
                           double& optimal_time);
  std::vector<double> cubic(double a, double b, double c, double d);
  std::vector<double> quadratic(double a, double b, double c);

  /* my one shot, use min acc integral */
  std::vector<double> quartic(double a, double b, double c, double d, double e);
  double getOptimalTime(Eigen::Vector3d p0, Eigen::Vector3d p1,
                        Eigen::Vector3d v0);
  Eigen::MatrixXd getShotTrajectory(Eigen::Vector3d p0, Eigen::Vector3d p1,
                                    Eigen::Vector3d v0, Eigen::Vector3d& v1,
                                    double& T);
  bool freeEndVelShot(GridNodePtr node1, GridNodePtr node2);

  /* searching */
  void stateTransit1(Eigen::VectorXd init_state, Eigen::VectorXd& final_state,
                     Eigen::Vector3d um, double tau);
  void stateTransitPrecompute(Eigen::VectorXd pre_state,
                              Eigen::VectorXd& final_state, Eigen::Vector3d um,
                              int step, double delta_t);
  void jerkInputStateTransit(Eigen::VectorXd init_state,
                             Eigen::VectorXd& final_state, Eigen::Vector3d um,
                             double tau);
  // diff, (edge_cost,heu_cost), (state, optimal_time)
  void getNeighbor(GridNodePtr current_node, GridNodePtr end_node,
                   Neighbors& neighbors, double& transit_time1,
                   double& transit_time2, double& vis_time);
  void getNeighborInit(GridNodePtr current_node, GridNodePtr end_node,
                       Neighbors& neighbors, double& transit_time1,
                       double& transit_time2, double& vis_time);

 public:
  HybridAStarPathFinder(ros::NodeHandle node) {
    // path publisher
    path_pub =
        node.advertise<visualization_msgs::Marker>("hybridastar/path", 10);
  };
  ~HybridAStarPathFinder(){};

  /* main API */
  void setParameterAuto();
  void resetNode();
  void setHorizon(double hrz) { this->horizon = hrz; }
  int searchPath(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel,
                 Eigen::Vector3d start_acc, Eigen::Vector3d end_pt,
                 Eigen::Vector3d end_vel, bool init);

  std::vector<GridNodePtr> retrievePath(GridNodePtr current);
  void initGridNodeMap(Eigen::Vector3i grid_num, double _resolution,
                       Eigen::Vector3d origin);
  void linkLocalMap(SDFMap& sdf_map);
  void resetLocalMap();
  void resetPath();
  std::vector<Eigen::Vector3d> getPath();
  std::vector<GridNodePtr> getVisitedNodes();
  std::vector<GridNodePtr> getPathNodes();
  std::vector<Eigen::Vector3d> getKinoTraj(double delta_t);
  Eigen::MatrixXd getSamples(double& ts, int& K, int N = 1,
                             bool repeat = false);
  void getKinoTrajMat(Eigen::MatrixXd& Pos, Eigen::MatrixXd& Vel,
                      Eigen::MatrixXd& Acc, Eigen::VectorXd& Time);

  // color 1234, rgby, resolution = 0.05
  ros::Publisher path_pub;
  void displayPathWithColor(std::vector<Eigen::Vector3d> path,
                            double resolution, int color, int id) {
    visualization_msgs::Marker mk;
    mk.header.frame_id = "world";
    mk.header.stamp = ros::Time::now();
    mk.type = visualization_msgs::Marker::SPHERE_LIST;
    mk.action = visualization_msgs::Marker::DELETE;
    mk.id = id;

    path_pub.publish(mk);

    mk.action = visualization_msgs::Marker::ADD;
    mk.pose.orientation.x = 0.0;
    mk.pose.orientation.y = 0.0;
    mk.pose.orientation.z = 0.0;
    mk.pose.orientation.w = 1.0;
    mk.color.a = 1.0;

    if (color == 1) {
      mk.color.r = 1.0;
      mk.color.g = 0.0;
      mk.color.b = 0.0;
    } else if (color == 2) {
      mk.color.g = 1.0;
      mk.color.r = 0.0;
      mk.color.b = 0.0;
    } else if (color == 3) {
      mk.color.b = 1.0;
      mk.color.g = 0.0;
      mk.color.r = 0.0;
    } else if (color == 4) {
      mk.color.g = 1.0;
      mk.color.r = 1.0;
      mk.color.b = 0.0;
    }

    mk.scale.x = resolution;
    mk.scale.y = resolution;
    mk.scale.z = resolution;

    geometry_msgs::Point pt;
    for (int i = 0; i < int(path.size()); i++) {
      pt.x = path[i](0);
      pt.y = path[i](1);
      pt.z = path[i](2);
      mk.points.push_back(pt);
    }
    path_pub.publish(mk);
    ros::Duration(0.001).sleep();
  }
};

#endif