#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>

#include <stdlib.h>
#include "grad_traj_optimization/display.h"
#include "grad_traj_optimization/grad_traj_optimizer.h"

using namespace std;

int main(int argc, char **argv) {
  ros::init(argc, argv, "random");
  ros::NodeHandle node;

  // -------------------ros initialization---------------------
  setpoint_pub =
      node.advertise<visualization_msgs::Marker>("trajopt/setpoint", 10);
  traj_point_pub =
      node.advertise<visualization_msgs::Marker>("trajopt/traj_point", 10);
  traj_pub = node.advertise<nav_msgs::Path>("trajopt/init_traj", 5);
  ros::Publisher visualization_pub = node.advertise<visualization_msgs::Marker>(
      "sdf_tools_tutorial_visualization", 1, true);

  srand(ros::Time::now().toSec());
  ros::Duration(0.5).sleep();

  // randomly generate some way points
  Eigen::Vector3d start, end;
  start(0) = start(1) = -4.0;
  end(0) = end(1) = 4.0;
  start(2) = end(2) = 2.0;

  vector<Eigen::Vector3d> way_points;
  way_points.push_back(start);

  double min_step = 1.5, max_step = 2.0;
  double dist = (start - end).norm();
  for (int i = 0; dist > min_step;) {
    Eigen::Vector3d wp;
    wp(0) = way_points[i](0) + max_step * rand() / double(RAND_MAX);
    wp(1) = way_points[i](1) + max_step * rand() / double(RAND_MAX);
    wp(2) = 2.0;

    if (fabs(wp(0)) <= 4.0 && fabs(wp(1)) <= 4.0 &&
        (wp - way_points[i]).norm() > min_step &&
        (wp - way_points[i]).norm() < max_step) {
      dist = (wp - end).norm();
      way_points.push_back(wp);
      ++i;
      // cout << "new way point:" << wp << endl;
    }
  }

  way_points.push_back(end);
  point_num = way_points.size();

  //  display  waypoints in rviz
  visualizeSetPoints(way_points);
  cout << "Way points created" << endl;

  int obs_num = 50;
  vector<Eigen::Vector3d> obstacles;
  cout << "Add obstacles to map" << endl;
  int fail_num = 0;
  for (int i = 0; i < obs_num;) {
    // randomly create a obstacle point
    Eigen::Vector3d pt;
    pt(0) = -3.5 + 7.0 * rand() / double(RAND_MAX);
    pt(1) = -3.5 + 7.0 * rand() / double(RAND_MAX);
    pt(2) = 2.0;

    // ensure that obstacle is far enough from start and end
    if ((pt - start).norm() < 2.0 || (pt - end).norm() < 2.0) continue;

    // ensure that obstacle is far enough from waypoint
    double dist_thresh = 1.1;
    double min_dist = 1000.0;
    for (int j = 0; j < way_points.size(); ++j) {
      double dist = (way_points[j] - pt).norm();
      if (dist < min_dist) min_dist = dist;
    }

    if (min_dist > dist_thresh) {
      // ensure that each obstacle is far enough from others
      if (i == 0) {
        obstacles.push_back(pt);
        ++i;
      } else {
        min_dist = 1000.0;
        dist_thresh = 1.2;
        for (int j = 0; j < obstacles.size(); ++j) {
          double dist = (obstacles[j] - pt).norm();
          if (dist < min_dist) min_dist = dist;
        }

        if (min_dist > dist_thresh) {
          obstacles.push_back(pt);
          ++i;
          fail_num = 0;
        } else {
          ++fail_num;
        }
      }
    }
    if (fail_num > 10000) {
      break;
    }
  }

  cout << "----------------------Obstacles generated!----------------------"
       << endl;

  // add the generated obstacles into collision map

  // ---------------------main optimization procedure---------------------
  GradTrajOptimizer grad_traj_opt(node, way_points);
  grad_traj_opt.initSDFMap(Eigen::Vector3d(40, 40, 5),
                           Eigen::Vector3d(-40 / 2, -40 / 2, 0.0), 0.2);
  grad_traj_opt.updateSDFMap(obstacles);

  Eigen::MatrixXd coeff;
  grad_traj_opt.getCoefficient(coeff);
  grad_traj_opt.getSegmentTime(my_time);
  displayTrajectory(coeff, false);

  // first step optimization
  grad_traj_opt.optimizeTrajectory(OPT_FIRST_STEP);
  grad_traj_opt.getCoefficient(coeff);
  displayTrajectory(coeff, false);

  // second step optimization
  grad_traj_opt.optimizeTrajectory(OPT_SECOND_STEP);
  grad_traj_opt.getCoefficient(coeff);
  displayTrajectory(coeff, true);

  return 0;
}
