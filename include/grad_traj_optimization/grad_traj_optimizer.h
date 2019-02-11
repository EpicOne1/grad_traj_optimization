#ifndef _GRAD_TRAJ_OPTIMIZER_H_
#define _GRAD_TRAJ_OPTIMIZER_H_

#include <Eigen/Eigen>
#include <nlopt.hpp>

// just use to get signed distance field
#include <ros/ros.h>
#include "grad_traj_optimization/sdf_map.h"

#include "qp_generator.h"

#define GDTB getDistanceToBoundary
#define OPT_INITIAL_TRY 0
#define OPT_FIRST_STEP 1
#define OPT_SECOND_STEP 2

using namespace std;

class GradTrajOptimizer {
 public:
  GradTrajOptimizer();

  void setPath(const vector<Eigen::Vector3d> &way_points);

  bool optimizeTrajectory(int step);

  void getCoefficient(Eigen::MatrixXd &coeff);

  void getSegmentTime(Eigen::VectorXd &seg_time);

  void initSDFMap(Eigen::Vector3d map_size_3d, Eigen::Vector3d origin,
                  double resolution);

  void updateSDFMap(vector<Eigen::Vector3d> obs);

 private:
  /** coefficient of polynomials*/
  Eigen::MatrixXd coeff;

  /** important matrix and  variables*/
  Eigen::MatrixXd A;
  Eigen::MatrixXd C;
  Eigen::MatrixXd L;
  Eigen::MatrixXd R;
  Eigen::MatrixXd Rff;
  Eigen::MatrixXd Rpp;
  Eigen::MatrixXd Rpf;
  Eigen::MatrixXd Rfp;
  Eigen::VectorXd segment_time;
  Eigen::MatrixXd V;
  Eigen::MatrixXd Df;
  Eigen::MatrixXd Dp;
  Eigen::MatrixXd origin_dp;
  Eigen::MatrixXd initial_dp;
  Eigen::MatrixXd path;  // nx3
  int num_dp;
  int num_df;
  int num_point;
  mutable int iter_num = 0;
  mutable double total_time = 0.0;
  int step = 1;
  double resolution;
  int algorithm;
  double time_limit_1, time_limit_2;
  double deltat;
  double bos, vos, aos;

  /** dynamics parameter from param server*/
  double w_smooth, w_collision;
  double d0, r, alpha, v0, r_v, alpha_v, a0, r_a, alpha_a;

  double mean_v, mean_a, init_time;

  /** optimizer*/
  nlopt::opt optimizer;

  /* sdf map */
  SDFMap sdf_map;

  /** main computation function,get smoothness, collision ,velocity ,accleration
   * cost and gradient*/
  void getCostAndGradient(std::vector<double> dp, double &cost,
                          std::vector<double> &grad);

  /** cost function of optimization */
  static double costFunc(const std::vector<double> &x,
                         std::vector<double> &grad, void *func_data);

  /** convert derivatives of end points to polynomials coefficient */
  void getCoefficientFromDerivative(Eigen::MatrixXd &coeff,
                                    const std::vector<double> &_dp);

  /** get distance and gradient in signed distance field ,by position query*/
  void getDistanceAndGradient(Eigen::Vector3d &pos, double &dist,
                              Eigen::Vector3d &grad);

  void getPositionFromCoeff(Eigen::Vector3d &pos, const Eigen::MatrixXd &coeff,
                            const int &index, const double &time);

  void getVelocityFromCoeff(Eigen::Vector3d &vel, const Eigen::MatrixXd &coeff,
                            const int &index, const double &time);

  void getAccelerationFromCoeff(Eigen::Vector3d &acc,
                                const Eigen::MatrixXd &coeff, const int &index,
                                const double &time);

  /** penalty and gradient */
  void getDistancePenalty(const double &distance, double &cost);
  void getDistancePenaltyGradient(const double &distance, double &grad);
  void getVelocityPenalty(const double &distance, double &cost);
  void getVelocityPenaltyGradient(const double &vel, double &grad);
  void getAccelerationPenalty(const double &distance, double &cost);
  void getAccelerationPenaltyGradient(const double &acc, double &grad);

  void getTimeMatrix(const double &t, Eigen::MatrixXd &T);
};

#endif