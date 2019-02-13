#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>
#include "std_msgs/Empty.h"

#include <stdlib.h>
// #include "grad_traj_optimization/a_star.h"
// #include <grad_traj_optimization/rrgPathFinder.h>
#include <grad_traj_optimization/path_finder.h>
#include "grad_traj_optimization/display.h"
#include "grad_traj_optimization/grad_traj_optimizer.h"

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include "sensor_msgs/PointCloud.h"

using namespace std;

pcl::PointCloud<pcl::PointXYZ> latest_cloud;

double start_x, start_y, start_z;
double start_vx, start_vy, start_vz;
double goal_x, goal_y, goal_z;

bool have_goal, have_map;

void mapCallback(const sensor_msgs::PointCloud2 &msg) {
  pcl::fromROSMsg(msg, latest_cloud);

  if ((int)latest_cloud.points.size() == 0) return;

  have_map = true;
  std::cout << "[1]: get map" << std::endl;
}

void sgCallback(const sensor_msgs::PointCloud &msg) {
  if (msg.points.size() != 3) {
    std::cout << "sg num error." << std::endl;
    return;
  } else {
    cout << "get sg msg." << endl;
  }

  start_x = msg.points[0].x, start_y = msg.points[0].y,
  start_z = msg.points[0].z;
  start_vx = msg.points[1].x, start_vy = msg.points[1].y,
  start_vz = msg.points[1].z;
  goal_x = msg.points[2].x, goal_y = msg.points[2].y, goal_z = msg.points[2].z;

  have_goal = true;
  std::cout << "[1]: get sg" << std::endl;
}

int main(int argc, char **argv) {
  ros::init(argc, argv, "random");
  ros::NodeHandle node;

  // -------------------ros initialization---------------------
  setpoint_pub =
      node.advertise<visualization_msgs::Marker>("grad_opt/setpoint", 10);
  traj_point_pub =
      node.advertise<visualization_msgs::Marker>("grad_opt/traj_point", 10);
  traj_pub = node.advertise<nav_msgs::Path>("grad_opt/init_traj", 5);

  ros::Subscriber map_sub =
      node.subscribe("/laser_cloud_surround", 1, mapCallback);
  ros::Subscriber sg_sub = node.subscribe("/start_goal", 1, sgCallback);
  ros::Publisher finish_pub =
      node.advertise<std_msgs::Empty>("/finish_test", 1, true);

  srand(ros::Time::now().toSec());
  ros::Duration(0.5).sleep();

  have_goal = false, have_map = false;
  const int use_map_num = 10;
  int exp_num = 0;

  /* init traj generator */
  GradTrajOptimizer grad_traj_opt;
  grad_traj_opt.initSDFMap(Eigen::Vector3d(40, 40, 5),
                           Eigen::Vector3d(-40 / 2, -40 / 2, 0.0), 0.2);

  /* init path finder */
  rrtPathFinder path_finder;
  path_finder.setParam(0.5, 0.2, 5.0, 10);

  /* main loop */
  while (ros::ok()) {
    /* wait for map and sg ready */
    while (ros::ok()) {
      if (have_map && have_goal) break;
      ros::Duration(0.1).sleep();
      ros::spinOnce();
    }
    cout << "[2]: Map and SG ok!" << endl;

    /* manage map */
    if (exp_num % use_map_num == 0) {
      vector<Eigen::Vector3d> obss;
      for (int i = 0; i < latest_cloud.points.size(); ++i) {
        pcl::PointXYZ pt = latest_cloud.points[i];
        obss.push_back(Eigen::Vector3d(pt.x, pt.y, pt.z));
      }
      grad_traj_opt.updateSDFMap(obss);
      path_finder.setInput(latest_cloud);
    }

    /* path finding */
    Eigen::Vector3d start, end;
    start(0) = start_x, start(1) = start_y, start(2) = start_z;
    end(0) = goal_x, end(1) = goal_y, end(2) = goal_z;

    Point point_s, point_g;
    point_s.x = start_x, point_s.y = start_y, point_s.z = start_z;
    point_g.x = goal_x, point_g.y = goal_y, point_g.z = goal_z;

    ros::Time t1 = ros::Time::now();

    path_finder.reset();
    path_finder.setPt(point_s, point_g, -20, 20, -20, 20, 0, 5, 0, 1, 10, 50000,
                      0.15, 0.05);
    path_finder.RRTpathFind(0.05);

    // vector<Eigen::Vector3d> path = path_finder.getPath();
    if (!path_finder.path_find_state) {
      ROS_WARN("[2]: can't find a path");
    } else {
      ros::Time t2 = ros::Time::now();
      double time_search = (t2 - t1).toSec();
      cout << "time in search: " << time_search << endl;

      Eigen::MatrixXd path_mat = path_finder.getPath().first;
      cout << "path:\n" << path_mat << endl;
      vector<Eigen::Vector3d> path;
      for (int i = 0; i < path_mat.rows(); ++i) {
        if (i == 1) continue;
        if (i == path_mat.rows() - 2) continue;
        Eigen::Vector3d wp = path_mat.row(i);
        path.push_back(wp);
      }

      if (path.size() == 2) {
        Eigen::Vector3d mid = 0.5 * (path[0] + path[1]);
        path.insert(path.begin() + 1, mid);
      }
      point_num = path.size();
      cout << "p num: " << point_num << endl;

      visualizeSetPoints(path);

      /* generate traj */
      t1 = ros::Time::now();
      grad_traj_opt.setPath(path);

      Eigen::MatrixXd coeff;
      grad_traj_opt.getSegmentTime(my_time);
      cout << "time: " << my_time.transpose() << endl;
      grad_traj_opt.getCoefficient(coeff);
      displayTrajectory(coeff, false);

      // second step optimization
      grad_traj_opt.optimizeTrajectory(OPT_SECOND_STEP);
      grad_traj_opt.getCoefficient(coeff);
      displayTrajectory(coeff, false);

      t2 = ros::Time::now();
      double time_opt = (t2 - t1).toSec();
      cout << "time in opt: " << time_opt << endl;

      cout << "total_time: " << time_search + time_opt << endl;
      cout << "[2]: finish test." << endl;
    }

    /* finish test flag */
    ++exp_num;
    have_goal = false;
    if (exp_num % use_map_num == 0) have_map = false;

    std_msgs::Empty finish_msg;
    finish_pub.publish(finish_msg);

    // ros::Duration(0.5).sleep();
  }

  return 0;
}
