#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>
#include "std_msgs/Empty.h"

#include <stdlib.h>
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
      node.advertise<visualization_msgs::Marker>("trajopt/setpoint", 10);
  traj_point_pub =
      node.advertise<visualization_msgs::Marker>("trajopt/traj_point", 10);
  traj_pub = node.advertise<nav_msgs::Path>("trajopt/init_traj", 5);
  ros::Publisher visualization_pub = node.advertise<visualization_msgs::Marker>(
      "sdf_tools_tutorial_visualization", 1, true);

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

  /* main loop */
  while (ros::ok()) {
    /* wait for map and sg ready */
    while (ros::ok()) {
      if (have_map && have_goal) break;
      ros::Duration(0.5).sleep();
      ros::spinOnce();
    }
    cout << "[1]: Map and SG ok!" << endl;

    // path finding
    Eigen::Vector3d start, end;
    start(0) = start_x, start(1) = start_y, start(2) = start_z;
    end(0) = goal_x, end(1) = goal_y, end(2) = goal_z;

    vector<Eigen::Vector3d> path;
    const double lambda = 0.2;
    for (int i = 0; i <= 5; ++i) {
      Eigen::Vector3d mid_pt = (1 - lambda * i) * start + lambda * i * end;
      path.push_back(mid_pt);
    }
    point_num = path.size();

    //  display  waypoints in rviz
    visualizeSetPoints(path);
    cout << "Way points created" << endl;

    /* convert map */
    vector<Eigen::Vector3d> obss;
    for (int i = 0; i < latest_cloud.points.size(); ++i) {
      pcl::PointXYZ pt = latest_cloud.points[i];
      obss.push_back(Eigen::Vector3d(pt.x, pt.y, pt.z));
    }

    // optimization procedure
    GradTrajOptimizer grad_traj_opt(node, path);
    grad_traj_opt.initSDFMap(Eigen::Vector3d(40, 40, 5),
                             Eigen::Vector3d(-40 / 2, -40 / 2, 0.0), 0.2);
    grad_traj_opt.updateSDFMap(obss);

    Eigen::MatrixXd coeff;
    grad_traj_opt.getCoefficient(coeff);
    grad_traj_opt.getSegmentTime(my_time);
    displayTrajectory(coeff, false);

    // first step optimization
    // grad_traj_opt.optimizeTrajectory(OPT_FIRST_STEP);
    // grad_traj_opt.getCoefficient(coeff);
    // displayTrajectory(coeff, false);

    // second step optimization
    grad_traj_opt.optimizeTrajectory(OPT_SECOND_STEP);
    grad_traj_opt.getCoefficient(coeff);
    displayTrajectory(coeff, false);

    /* finish test flag */
    ++exp_num;
    have_goal = false;
    if (exp_num % use_map_num == 0) have_map = false;

    std_msgs::Empty finish_msg;
    finish_pub.publish(finish_msg);
  }

  return 0;
}
