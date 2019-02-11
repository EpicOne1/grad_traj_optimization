#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>
#include "std_msgs/Empty.h"

#include <stdlib.h>
#include "grad_traj_optimization/a_star.h"
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

void simplifyPath(vector<Eigen::Vector3d> &path, SDFMap &sdf_map) {
  vector<Eigen::Vector3d> path_sp;
  path_sp.push_back(path.front());

  int cur_id = 1;
  while (cur_id < path.size()) {
    Eigen::Vector3d pt1 = path_sp.back();
    Eigen::Vector3d pt2 = path[cur_id];

    /* detect any collision between pt1 and pt2 */
    bool collision = false;
    Eigen::Vector3d dir = pt2 - pt1;
    double lambda = 0.1;
    while (lambda < 1.0) {
      Eigen::Vector3d ptm = pt1 + lambda * dir;
      // int occ = sdf_map.getOccupancy(ptm);
      double dist = sdf_map.getDistance(ptm);
      if (dist < 0.4) {
        collision = true;
        break;
      }
      lambda += 0.1;
    }

    if (!collision) {
      ++cur_id;
    } else {
      path_sp.push_back(path[cur_id - 1]);
    }
  }
  path_sp.push_back(path.back());

  path = path_sp;
}

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

  /* init traj generator */
  GradTrajOptimizer grad_traj_opt;
  grad_traj_opt.initSDFMap(Eigen::Vector3d(40, 40, 5),
                           Eigen::Vector3d(-40 / 2, -40 / 2, 0.0), 0.2);

  /* init path finder */
  AStarPathFinder *path_finder = new AStarPathFinder();
  Eigen::Vector3i gn;
  gn(0) = ceil(40 / 0.2), gn(1) = ceil(40 / 0.2), gn(2) = ceil(5 / 0.2);
  path_finder->initGridNodeMap(gn, 0.2, Eigen::Vector3d(-40 / 2, -40 / 2, 0.0));

  /* init sdf */
  SDFMap sdf_map(Eigen::Vector3d(-40 / 2, -40 / 2, 0.0), 0.2,
                 Eigen::Vector3d(40, 40, 5));

  /* main loop */
  while (ros::ok()) {
    /* wait for map and sg ready */
    while (ros::ok()) {
      if (have_map && have_goal) break;
      ros::Duration(0.1).sleep();
      ros::spinOnce();
    }
    cout << "[1]: Map and SG ok!" << endl;

    /* manage map */
    if (exp_num % use_map_num == 0) {
      vector<Eigen::Vector3d> obss;
      for (int i = 0; i < latest_cloud.points.size(); ++i) {
        pcl::PointXYZ pt = latest_cloud.points[i];
        obss.push_back(Eigen::Vector3d(pt.x, pt.y, pt.z));
      }
      grad_traj_opt.updateSDFMap(obss);

      sdf_map.resetBuffer();
      for (int i = 0; i < int(obss.size()); ++i) {
        sdf_map.setOccupancy(obss[i]);
      }
      sdf_map.updateESDF3d();

      path_finder->linkLocalMap(sdf_map);
    }

    // path finding
    Eigen::Vector3d start, end;
    start(0) = start_x, start(1) = start_y, start(2) = start_z;
    end(0) = goal_x, end(1) = goal_y, end(2) = goal_z;

    path_finder->resetNode();
    path_finder->resetPath();
    path_finder->searchPath(start, end);

    vector<Eigen::Vector3d> path = path_finder->getPath();

    /* simplify path */
    cout << "path size bf sp: " << path.size() << endl;
    simplifyPath(path, sdf_map);
    if (path.size() == 2) {
      path.clear();
      const double lambda = 0.2;
      for (int i = 0; i <= 5; ++i) {
        Eigen::Vector3d mid_pt = (1 - lambda * i) * start + lambda * i * end;
        path.push_back(mid_pt);
      }
    }
    cout << "path size af sp: " << path.size() << endl;

    point_num = path.size();

    /* straight line initialization */
    // const double lambda = 0.2;
    // for (int i = 0; i <= 5; ++i) {
    //   Eigen::Vector3d mid_pt = (1 - lambda * i) * start + lambda * i * end;
    //   path.push_back(mid_pt);
    // }

    //  display  waypoints in rviz
    visualizeSetPoints(path);
    cout << "Way points created" << endl;

    /* generate traj */
    grad_traj_opt.setPath(path);

    Eigen::MatrixXd coeff;
    // grad_traj_opt.getCoefficient(coeff);
    // displayTrajectory(coeff, false);

    // // first step optimization
    // grad_traj_opt.optimizeTrajectory(OPT_FIRST_STEP);
    // grad_traj_opt.getCoefficient(coeff);
    // displayTrajectory(coeff, false);

    // second step optimization
    grad_traj_opt.optimizeTrajectory(OPT_SECOND_STEP);
    grad_traj_opt.getCoefficient(coeff);
    grad_traj_opt.getSegmentTime(my_time);
    displayTrajectory(coeff, false);

    /* finish test flag */
    ++exp_num;
    have_goal = false;
    if (exp_num % use_map_num == 0) have_map = false;

    std_msgs::Empty finish_msg;
    finish_pub.publish(finish_msg);

    ros::Duration(0.5).sleep();
  }

  return 0;
}
