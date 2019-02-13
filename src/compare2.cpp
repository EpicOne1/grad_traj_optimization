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

ros::Publisher poly_traj_pub;

void displayPathWithColor(vector<Eigen::Vector3d> path, double resolution,
                          Eigen::Vector4d color, int id) {
  visualization_msgs::Marker mk;
  mk.header.frame_id = "world";
  mk.header.stamp = ros::Time::now();
  mk.type = visualization_msgs::Marker::SPHERE_LIST;
  mk.action = visualization_msgs::Marker::DELETE;
  mk.id = id;
  poly_traj_pub.publish(mk);

  mk.action = visualization_msgs::Marker::ADD;
  mk.pose.orientation.x = 0.0, mk.pose.orientation.y = 0.0,
  mk.pose.orientation.z = 0.0, mk.pose.orientation.w = 1.0;
  mk.scale.x = resolution, mk.scale.y = resolution, mk.scale.z = resolution;

  mk.color.r = color(0), mk.color.g = color(1), mk.color.b = color(2),
  mk.color.a = color(3);

  geometry_msgs::Point pt;
  for (int i = 0; i < int(path.size()); i++) {
    pt.x = path[i](0), pt.y = path[i](1), pt.z = path[i](2);
    mk.points.push_back(pt);
  }
  poly_traj_pub.publish(mk);
  ros::Duration(0.01).sleep();
}

class PolynomialTraj {
 private:
  vector<double> times;        // time of each segment
  vector<vector<double>> cxs;  // coefficient of x of each segment
  vector<vector<double>> cys;  // coefficient of y of each segment
  vector<vector<double>> czs;  // coefficient of z of each segment

  double time_sum;
  int num_seg;

 public:
  PolynomialTraj(/* args */) {}
  ~PolynomialTraj() {}

  void reset() {
    times.clear(), cxs.clear(), cys.clear(), czs.clear();
    time_sum = 0.0, num_seg = 0;
  }

  void addSegment(vector<double> cx, vector<double> cy, vector<double> cz,
                  double t) {
    cxs.push_back(cx), cys.push_back(cy), czs.push_back(cz), times.push_back(t);
  }

  void init() {
    num_seg = times.size();
    time_sum = 0.0;
    for (int i = 0; i < times.size(); ++i) {
      time_sum += times[i];
    }
  }

  double getTimeSum() { return this->time_sum; }

  Eigen::Vector3d evaluate(double t) {
    /* detetrmine segment num */
    int idx = 0;
    while (times[idx] < t) {
      t -= times[idx];
      ++idx;
    }

    /* evaluation */
    int order = cxs[idx].size();
    Eigen::VectorXd cx(order), cy(order), cz(order), tv(order);
    for (int i = 0; i < order; ++i) {
      cx(i) = cxs[idx][i], cy(i) = cys[idx][i], cz(i) = czs[idx][i];
      tv(order - 1 - i) = std::pow(t, double(i));
    }

    Eigen::Vector3d pt;
    pt(0) = tv.dot(cx), pt(1) = tv.dot(cy), pt(2) = tv.dot(cz);
    return pt;
  }
};

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
  // traj_pub = node.advertise<nav_msgs::Path>("grad_opt/init_traj", 5);
  poly_traj_pub = node.advertise<visualization_msgs::Marker>(
      "/gradient_based/traj", 10, true);

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
      Eigen::VectorXd time_sgm;
      grad_traj_opt.getSegmentTime(time_sgm);
      my_time = time_sgm;
      // grad_traj_opt.getCoefficient(coeff);
      // displayTrajectory(coeff, false);

      // second step optimization
      grad_traj_opt.optimizeTrajectory(OPT_SECOND_STEP);
      grad_traj_opt.getCoefficient(coeff);

      t2 = ros::Time::now();
      double time_opt = (t2 - t1).toSec();
      cout << "time in opt: " << time_opt << endl;
      cout << "total_time: " << time_search + time_opt << endl;

      /* convert coefficient to poly_traj */
      PolynomialTraj poly_traj;
      for (int i = 0; i < coeff.rows(); ++i) {
        vector<double> cx(6), cy(6), cz(6);
        for (int j = 0; j < 6; ++j) {
          cx[j] = coeff(i, j), cy[j] = coeff(i, j + 6),
          cz[j] = coeff(i, j + 12);
        }
        reverse(cx.begin(), cx.end());
        reverse(cy.begin(), cy.end());
        reverse(cz.begin(), cz.end());
        double ts = time_sgm(i);
        poly_traj.addSegment(cx, cy, cz, ts);
      }
      poly_traj.init();

      double time_sum = poly_traj.getTimeSum();
      double eval_t = 0.0;
      vector<Eigen::Vector3d> traj_vis;
      while (eval_t < time_sum) {
        Eigen::Vector3d pt = poly_traj.evaluate(eval_t);
        traj_vis.push_back(pt);
        eval_t += 0.01;
      }
      displayPathWithColor(traj_vis, 0.1, Eigen::Vector4d(0, 1, 0, 1), 1);

      // displayTrajectory(coeff, false);

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
