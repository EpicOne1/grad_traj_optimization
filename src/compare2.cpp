#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>
#include "std_msgs/Empty.h"

#include <stdlib.h>
// #include "grad_traj_optimization/a_star.h"
// #include <grad_traj_optimization/rrgPathFinder.h>
#include <grad_traj_optimization/hybrid_astar.h>
#include "grad_traj_optimization/display.h"
#include "grad_traj_optimization/grad_traj_optimizer.h"

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <fstream>
#include <grad_traj_optimization/polynomial_traj.hpp>
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
  unique_ptr<HybridAStarPathFinder> path_finder;
  path_finder.reset(new HybridAStarPathFinder(node));
  path_finder->setParameterAuto();  // need to modify launch

  Eigen::Vector3i gn;
  gn(0) = ceil(40 / 0.2), gn(1) = ceil(40 / 0.2), gn(2) = ceil(5.0 / 0.2);
  path_finder->initGridNodeMap(gn, 0.2, Eigen::Vector3d(-40 / 2, -40 / 2, 0.0));

  SDFMap sdf_map(Eigen::Vector3d(-40 / 2, -40 / 2, 0.0), 0.2,
                 Eigen::Vector3d(40, 40, 5));

  /* main loop */
  while (ros::ok()) {
    /* wait for map and sg ready */
    while (ros::ok()) {
      if (have_map && have_goal) break;
      // ros::Duration(0.1).sleep();
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

      sdf_map.resetBuffer(Eigen::Vector3d(-40 / 2, -40 / 2, 0.0),
                          Eigen::Vector3d(40 / 2, 40 / 2, 5.0));
      for (int i = 0; i < int(obss.size()); ++i) {
        sdf_map.setOccupancy(obss[i]);
      }
      sdf_map.updateESDF3d();
      path_finder->linkLocalMap(sdf_map);
    }

    /* manage start and goal */
    Eigen::Vector3d start_pt, start_vel, start_acc, end_pt, end_vel;
    start_pt(0) = start_x, start_pt(1) = start_y, start_pt(2) = start_z;
    // start_vel(0) = start_vx, start_vel(1) = start_vy, start_vel(2) =
    // start_vz;
    end_pt(0) = goal_x, end_pt(1) = goal_y, end_pt(2) = goal_z;
    start_vel.setZero();
    start_acc.setZero(), end_vel.setZero();

    /* path finding */
    ros::Time t1, t2;

    t1 = ros::Time::now();
    path_finder->resetNode();
    path_finder->resetPath();

    int status = path_finder->searchPath(start_pt, start_vel, start_acc, end_pt,
                                         end_vel, false);

    t2 = ros::Time::now();
    double time_search = (t2 - t1).toSec();
    cout << "time in search: " << time_search << endl;

    if (status == NO_PATH) {
      cout << "[2]:no path" << endl;
    } else {
      Eigen::MatrixXd Pos, Vel, Acc;
      Eigen::VectorXd Time;
      path_finder->getKinoTrajMat(Pos, Vel, Acc, Time);

      // visualizeSetPoints(path);

      /* ============================== front end evaluation ============================== */
      PolynomialTraj poly_traj1;
      double acc_cost2 = 0.0;
      for (int i = 0; i < Time.rows(); ++i) {
        vector<double> cx(6), cy(6), cz(6);
        std::fill(cx.begin(), cx.begin() + 3, 0.0);
        std::fill(cy.begin(), cy.begin() + 3, 0.0);
        std::fill(cz.begin(), cz.begin() + 3, 0.0);
        cx[3] = Acc(i, 0) / 2.0, cy[3] = Acc(i, 1) / 2.0,
        cz[3] = Acc(i, 2) / 2.0;
        cx[4] = Vel(i, 0), cy[4] = Vel(i, 1), cz[4] = Vel(i, 2);
        cx[5] = Pos(i, 0), cy[5] = Pos(i, 1), cz[5] = Pos(i, 2);

        poly_traj1.addSegment(cx, cy, cz, Time(i));

        /* acc cost */
        acc_cost2 += Acc.row(i).squaredNorm() * Time(i);
      }
      poly_traj1.init();
      cout << "[2]: acc cost for compare: " << acc_cost2 << endl;

      vector<Eigen::Vector3d> traj_vis_path = poly_traj1.getTraj();
      displayPathWithColor(traj_vis_path, 0.1, Eigen::Vector4d(1, 1, 0, 1),
      2);

      double acc_cost = poly_traj1.getAccCost();
      double time_sum_path = poly_traj1.getTimeSum();
      cout << "test2:" << exp_num + 1 << "solve_time:" << time_search
           << ",traj_time:" << time_sum_path << ",acc_cost:" << acc_cost
           << endl;

      std::ofstream file(
          "/home/bzhouai/workspaces/plan_ws/src/uav_planning_bm/resource/"
          "front2.txt",
          std::ios::app);
      if (file.fail()) {
        cout << "open file error!\n";
        return -1;
      }
      file << "test2:" << exp_num + 1 << "solve_time:" << time_search
           << ",traj_time:" << time_sum_path << ",acc_cost:" << acc_cost
           << "\n";
      file.close();

    /* ============================== compare optimization ============================== */
      bool compare_optimize = false;
      if (compare_optimize) {
        t1 = ros::Time::now();
        grad_traj_opt.setKinoPath(Pos, Vel, Acc, Time);
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
        vector<Eigen::Vector3d> traj_vis = poly_traj.getTraj();
        displayPathWithColor(traj_vis, 0.2, Eigen::Vector4d(0, 0, 1, 1), 1);

        /* evaluation */
        double time_sum, length, mean_v, max_v, mean_a, max_a, jerk;
        time_sum = poly_traj.getTimeSum();
        length = poly_traj.getLength();
        jerk = poly_traj.getJerk();
        poly_traj.getMeanAndMaxVel(mean_v, max_v);
        poly_traj.getMeanAndMaxAcc(mean_a, max_a);

        vector<double> vec_time, vec_cost;
        grad_traj_opt.getCostCurve(vec_cost, vec_time);

        cout << "test2:" << exp_num + 1 << ",jerk:" << jerk;
        cout << ",time:";
        for (int i = 0; i < vec_time.size(); ++i) {
          cout << vec_time[i];
          if (i != vec_time.size() - 1) cout << ";";
        }
        cout << ",cost:";
        for (int i = 0; i < vec_cost.size(); ++i) {
          cout << vec_cost[i];
          if (i != vec_cost.size() - 1)
            cout << ";";
          else
            cout << "\n";
        }
        std::ofstream file(
            "/home/bzhouai/workspaces/plan_ws/src/uav_planning_bm/resource/"
            "back2.txt",
            std::ios::app);
        if (file.fail()) {
          cout << "open file error!\n";
          return -1;
        }

        file << "test2:" << exp_num + 1 << ",jerk:" << jerk;
        file << ",time:";
        for (int i = 0; i < vec_time.size(); ++i) {
          file << vec_time[i];
          if (i != vec_time.size() - 1) file << ";";
        }
        file << ",cost:";
        for (int i = 0; i < vec_cost.size(); ++i) {
          file << vec_cost[i];
          if (i != vec_cost.size() - 1)
            file << ";";
          else
            file << "\n";
        }

        file.close();
      } /* optimize traj */

      // displayTrajectory(coeff, false);
    }
    /* ============================== compare optimization ============================== */

    /* finish test flag */
    cout << "[2]: finish test." << endl;
    ++exp_num;
    have_goal = false;
    if (exp_num % use_map_num == 0) have_map = false;

    ros::Duration(0.5).sleep();

    std_msgs::Empty finish_msg;
    finish_pub.publish(finish_msg);

    // ros::Duration(0.5).sleep();
  }

  return 0;
}
