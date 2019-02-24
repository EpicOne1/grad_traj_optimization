#include "grad_traj_optimization/grad_traj_optimizer.h"

GradTrajOptimizer::GradTrajOptimizer() {
  //-------------------------get parameter from server--------------------
  ros::param::get("/traj_opti_node1/alg", this->algorithm);
  ros::param::get("/traj_opti_node1/time_limit_1", this->time_limit_1);
  ros::param::get("/traj_opti_node1/time_limit_2", this->time_limit_2);
  ros::param::get("/traj_opti_node1/dt", this->deltat);

  ros::param::get("/traj_opti_node1/ws", this->w_smooth);
  ros::param::get("/traj_opti_node1/wc", this->w_collision);

  ros::param::get("/traj_opti_node1/alpha", this->alpha);
  ros::param::get("/traj_opti_node1/r", this->r);
  ros::param::get("/traj_opti_node1/d0", this->d0);

  ros::param::get("/traj_opti_node1/alpha_v", this->alpha_v);
  ros::param::get("/traj_opti_node1/r_v", this->r_v);
  ros::param::get("/traj_opti_node1/v0", this->v0);

  ros::param::get("/traj_opti_node1/alpha_a", this->alpha_a);
  ros::param::get("/traj_opti_node1/r_a", this->r_a);
  ros::param::get("/traj_opti_node1/a0", this->a0);

  /* upper and lower bound for optimization */
  ros::param::get("/traj_opti_node1/bos", this->bos);
  ros::param::get("/traj_opti_node1/vos", this->vos);
  ros::param::get("/traj_opti_node1/aos", this->aos);

  ros::param::get("/traj_opti_node1/mean_v", mean_v);
  ros::param::get("/traj_opti_node1/mean_a", mean_a);
  ros::param::get("/traj_opti_node1/init_time", init_time);
}

void GradTrajOptimizer::setKinoPath(Eigen::MatrixXd &Pos, Eigen::MatrixXd &Vel,
                                    Eigen::MatrixXd &Acc,
                                    Eigen::VectorXd &Time) {
  path = Pos;

  segment_time = Time;

  TrajectoryGenerator generator;
  // coeff = generator.PolyQPGeneration(path, vel, acc, segment_time, type);
  coeff = generator.PolyKinoGeneration(Pos, Vel, Acc, Time, 1);
  generator.StackOptiDep();
  R = generator.getR();
  Rff = generator.getRff();
  Rpp = generator.getRpp();
  Rpf = generator.getRpf();
  Rfp = generator.getRfp();
  L = generator.getL();
  A = generator.getA();
  C = generator.getC();

  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> d = generator.getInitialD();
  initial_dp = origin_dp = Dp = d.first;
  Df = d.second;

  V.resize(6, 6);
  for (int i = 0; i < 5; ++i) V(i, i + 1) = i + 1;

  num_dp = Dp.cols();
  num_df = Df.cols();
  num_point = segment_time.rows() + 1;
}

void GradTrajOptimizer::setPath(const vector<Eigen::Vector3d> &way_points) {
  /* generate optimization dependency */
  path = Eigen::MatrixXd::Zero(way_points.size(), 3);
  for (int i = 0; i < way_points.size(); ++i)
    path.row(i) = way_points[i].transpose();

  segment_time = Eigen::VectorXd::Zero(way_points.size() - 1);
  for (int i = 0; i < segment_time.size(); ++i) {
    double len = (path.row(i) - path.row(i + 1)).norm();
    if (i == 0 || i == segment_time.size()) {
      segment_time(i) = len / mean_v + init_time;
    } else {
      segment_time(i) = len / mean_v;
    }
  }

  Eigen::Vector3d vel, acc;
  vel.setZero();
  acc.setZero();
  int type = 2;

  TrajectoryGenerator generator;
  coeff = generator.PolyQPGeneration(path, vel, acc, segment_time, type);
  generator.StackOptiDep();
  R = generator.getR();
  Rff = generator.getRff();
  Rpp = generator.getRpp();
  Rpf = generator.getRpf();
  Rfp = generator.getRfp();
  L = generator.getL();
  A = generator.getA();
  C = generator.getC();

  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> d = generator.getInitialD();
  initial_dp = origin_dp = Dp = d.first;
  Df = d.second;

  V.resize(6, 6);
  for (int i = 0; i < 5; ++i) V(i, i + 1) = i + 1;

  num_dp = Dp.cols();
  num_df = Df.cols();
  num_point = segment_time.rows() + 1;
}

void GradTrajOptimizer::initSDFMap(Eigen::Vector3d map_size_3d,
                                   Eigen::Vector3d origin, double resolution) {
  sdf_map = SDFMap(origin, resolution, map_size_3d);
}

void GradTrajOptimizer::updateSDFMap(vector<Eigen::Vector3d> obs) {
  sdf_map.resetBuffer();

  for (int i = 0; i < int(obs.size()); ++i) {
    sdf_map.setOccupancy(obs[i]);
  }

  sdf_map.updateESDF3d();
  /* TODO: add virtual boundary */
}

bool GradTrajOptimizer::optimizeTrajectory(int step) {
  if (step != 0 && step != 1 && step != 2) {
    cout << "step number error, step should be 0, 1 or 2" << endl;
  }
  this->step = step;

  /* initilize nlopt */
  int seed = ros::Time::now().toNSec() % 65536;
  nlopt::srand(seed);
  nlopt::opt opt(nlopt::algorithm(this->algorithm),
                 3 * num_dp);  // x,y,z (3*n-3) x 3
  optimizer = opt;
  optimizer.set_min_objective(GradTrajOptimizer::costFunc, this);
  // optimizer.set_xtol_abs(1e-7);

  /* step specific option */
  if (step == OPT_FIRST_STEP) {
    optimizer.set_maxtime(time_limit_1);
  } else if (step == OPT_SECOND_STEP) {
    optimizer.set_maxtime(time_limit_2);
  }

  /* upper and lower bound */
  vector<double> lb, ub;
  lb.resize(3 * num_dp);
  ub.resize(3 * num_dp);
  for (int i = 0; i < num_dp; ++i) {
    if (i % 3 == 0) {
      lb[i] = path(i / 3 + 1, 0) - bos;
      lb[i + num_dp] = path(i / 3 + 1, 1) - bos;
      lb[i + num_dp * 2] = path(i / 3 + 1, 2) - bos;
      ub[i] = path(i / 3 + 1, 0) + bos;
      ub[i + num_dp] = path(i / 3 + 1, 1) + bos;
      ub[i + num_dp * 2] = path(i / 3 + 1, 2) + bos;
    } else if (i % 3 == 1) {
      lb[i] = -vos;
      lb[i + num_dp] = -vos;
      lb[i + 2 * num_dp] = -vos;
      ub[i] = vos;
      ub[i + num_dp] = vos;
      ub[i + num_dp * 2] = vos;
    } else {
      lb[i] = -aos;
      lb[i + num_dp] = -aos;
      lb[i + 2 * num_dp] = -aos;
      ub[i] = aos;
      ub[i + num_dp] = aos;
      ub[i + num_dp * 2] = aos;
    }
  }
  optimizer.set_lower_bounds(lb);
  optimizer.set_upper_bounds(ub);

  /* set initial value */
  std::vector<double> _dp(3 * num_dp);
  for (int i = 0; i < num_dp; ++i) {
    _dp[i] = Dp(0, i);
    _dp[i + num_dp] = Dp(1, i);
    _dp[i + 2 * num_dp] = Dp(2, i);
  }
  double min_f;

  /* optimize */
  cout << "-------------------begin optimization-------------------" << endl;
  vec_time.clear();
  vec_cost.clear();
  time_start = ros::Time::now();
  nlopt::result result = optimizer.optimize(_dp, min_f);

  // ---------------------------display the result---------------------------
  cout << "Optimized result is:" << result << endl;

  // ---------------------------update optimized
  // derivative---------------------------
  Dp.setZero();
  for (int i = 0; i < num_dp; ++i) {
    Dp(0, i) = _dp[i];
    Dp(1, i) = _dp[i + num_dp];
    Dp(2, i) = _dp[i + 2 * num_dp];
  }

  /* reallocate segment time */
  for (int i = 0; i < segment_time.size(); ++i) {
    double len = 0.0;
    // head and tail segment length
    if (i == 0) {
      len = sqrt(pow(Df(0, 0) - Dp(0, 0), 2) + pow(Df(1, 0) - Dp(1, 0), 2) +
                 pow(Df(2, 0) - Dp(2, 0), 2));
    } else if (i == segment_time.size() - 1) {
      len = sqrt(pow(Df(0, 3) - Dp(0, 3 * (i - 1)), 2) +
                 pow(Df(1, 3) - Dp(1, 3 * (i - 1)), 2) +
                 pow(Df(2, 3) - Dp(2, 3 * (i - 1)), 2));
    } else {
      len = sqrt(pow(Dp(0, 3 * (i - 1)) - Dp(0, 3 * i), 2) +
                 pow(Dp(1, 3 * (i - 1)) - Dp(1, 3 * i), 2) +
                 pow(Dp(2, 3 * (i - 1)) - Dp(2, 3 * i), 2));
    }
    segment_time(i) = len / mean_v;
    if (i == 0 || i == segment_time.size() - 1) segment_time(i) += init_time;
  }

  /* update optimized coefficient */
  getCoefficientFromDerivative(this->coeff, _dp);

  /* show optimization time */
  if (step == 1 || step == 2) {
    cout << "total time:" << total_time << endl
         << "iterative num:" << iter_num << endl;
    if (step == 2) {
      total_time = 0;
      iter_num = 0;
    }
  }

  return true;
}

void GradTrajOptimizer::getCoefficient(Eigen::MatrixXd &coe) {
  coe = this->coeff;
}

void GradTrajOptimizer::getSegmentTime(Eigen::VectorXd &seg_time) {
  seg_time = this->segment_time;
}

void GradTrajOptimizer::getCoefficientFromDerivative(
    Eigen::MatrixXd &coefficient, const std::vector<double> &_dp) {
  coefficient.resize(num_point - 1, 18);

  for (int i = 0; i < 3; ++i) {
    // merge df and dp ->d(df,dp)
    Eigen::VectorXd df(num_df);
    Eigen::VectorXd dp(num_dp);
    Eigen::VectorXd d(num_df + num_dp);

    df = Df.row(i);
    for (int j = 0; j < num_dp; j++) {
      dp(j) = _dp[j + num_dp * i];
    }

    d.segment(0, 6) = df;
    d.segment(6, num_dp) = dp;

    // convert derivative to coefficient
    Eigen::VectorXd coe(6 * (num_point - 1));
    coe = L * d;

    for (int j = 0; j < (num_point - 1); j++) {
      coefficient.block(j, 6 * i, 1, 6) = coe.segment(6 * j, 6).transpose();
    }
  }
}

void GradTrajOptimizer::getCostAndGradient(std::vector<double> dp, double &cost,
                                           std::vector<double> &_grad) {
  // get total iterative number and time
  iter_num++;
  ros::Time tb1 = ros::Time::now();

  // initialize
  cost = 0;
  double cost_smooth = 0;
  double cost_colli = 0;
  double cost_vel = 0;
  double cost_acc = 0;

  Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(3, num_dp);
  Eigen::MatrixXd g_smooth = Eigen::MatrixXd::Zero(3, num_dp);
  Eigen::MatrixXd g_colli = Eigen::MatrixXd::Zero(3, num_dp);
  Eigen::MatrixXd g_vel = Eigen::MatrixXd::Zero(3, num_dp);
  Eigen::MatrixXd g_acc = Eigen::MatrixXd::Zero(3, num_dp);

  // get smoothness cost
  // merge df and dp into d(df,dp
  Eigen::VectorXd dfx = Df.block(0, 0, 1, 6).transpose();
  Eigen::VectorXd dfy = Df.block(1, 0, 1, 6).transpose();
  Eigen::VectorXd dfz = Df.block(2, 0, 1, 6).transpose();

  Eigen::VectorXd dpx = Eigen::VectorXd::Zero(num_dp);
  Eigen::VectorXd dpy = Eigen::VectorXd::Zero(num_dp);
  Eigen::VectorXd dpz = Eigen::VectorXd::Zero(num_dp);
  for (int i = 0; i < num_dp; ++i) {
    dpx(i) = dp[i];
    dpy(i) = dp[i + num_dp];
    dpz(i) = dp[i + 2 * num_dp];
  }

  Eigen::VectorXd dx = Eigen::VectorXd::Zero(num_dp + num_df);
  Eigen::VectorXd dy = Eigen::VectorXd::Zero(num_dp + num_df);
  Eigen::VectorXd dz = Eigen::VectorXd::Zero(num_dp + num_df);
  dx.segment(0, 6) = dfx;
  dx.segment(6, num_dp) = dpx;
  dy.segment(0, 6) = dfy;
  dy.segment(6, num_dp) = dpy;
  dz.segment(0, 6) = dfz;
  dz.segment(6, num_dp) = dpz;

  // get smoothness cost,fs = d'Rd
  cost_smooth = double(dx.transpose() * R * dx) +
                double(dy.transpose() * R * dy) + (dz.transpose() * R * dz);

  // get smoothness gradient
  Eigen::MatrixXd gx_smooth = 2 * Rfp.transpose() * dfx + 2 * Rpp * dpx;
  Eigen::MatrixXd gy_smooth = 2 * Rfp.transpose() * dfy + 2 * Rpp * dpy;
  Eigen::MatrixXd gz_smooth = 2 * Rfp.transpose() * dfz + 2 * Rpp * dpz;

  g_smooth.row(0) = gx_smooth.transpose();
  g_smooth.row(1) = gy_smooth.transpose();
  g_smooth.row(2) = gz_smooth.transpose();

  // get polynomials coefficient, for evaluating penalty
  Eigen::MatrixXd coe;
  getCoefficientFromDerivative(coe, dp);

  // get coolision, velocity and acceleration cost and gradient
  Eigen::MatrixXd Ldp(6, num_dp);

  for (int s = 0; s < segment_time.size(); s++) {
    if (fabs(w_collision) < 1e-4) break;
    // get matrix Ldp
    Ldp = L.block(6 * s, 6, 6, num_dp);

    // discrete time step
    double dt = segment_time(s) / 30.0;

    for (double t = 1e-3; t < segment_time(s); t += dt) {
      // get position,velocity
      Eigen::Vector3d pos, vel;
      getPositionFromCoeff(pos, coe, s, t);
      getVelocityFromCoeff(vel, coe, s, t);
      double vel_norm = vel.norm() + 1e-5;

      // get information from sdf
      double dist = 0, gd = 0, cd = 0;
      Eigen::Vector3d grad;
      getDistanceAndGradient(pos, dist, grad);

      getDistancePenalty(dist, cd);
      getDistancePenaltyGradient(dist, gd);

      // time Matrix T
      Eigen::MatrixXd T(1, 6);
      getTimeMatrix(t, T);

      // collision cost
      cost_colli += cd * vel_norm * dt;

      // gradient of collision
      for (int k = 0; k < 3; k++) {
        g_colli.row(k) =
            g_colli.row(k) + (gd * grad(k) * cd * vel_norm * T * Ldp +
                              cd * (vel(k) / vel_norm) * T * V * Ldp) *
                                 dt;
      }
      // get velocity and accleration cost
      // if (step == 2) {
      //   double cv = 0, ca = 0, gv = 0, ga = 0;
      //   Eigen::Vector3d acc;
      //   getAccelerationFromCoeff(acc, coe, s, t);
      //   for (int k = 0; k < 3; k++) {
      //     getVelocityPenalty(vel(k), cv);
      //     cost_vel += cv * vel_norm * dt;

      //     getAccelerationPenalty(acc(k), ca);
      //     cost_acc += ca * vel_norm * dt;
      //   }
      //   for (int k = 0; k < 3; k++) {
      //     getVelocityPenaltyGradient(vel(k), gv);
      //     g_vel.row(k) =
      //         g_vel.row(k) + (gv * vel_norm * T * V * Ldp +
      //                         cv * (vel(k) / vel_norm) * T * V * Ldp) *
      //                            dt;

      //     getAccelerationPenaltyGradient(acc(k), ga);
      //     g_acc.row(k) =
      //         g_acc.row(k) + (ga * vel_norm * T * V * V * Ldp +
      //                         ca * (vel(k) / vel_norm) * T * V * Ldp) *
      //                            dt;
      //   }
      // }
    }
  }

  // sum up all cost
  double ws = this->w_smooth, wc = this->w_collision, wv = 1.0, wa = 1.0;
  if (step == OPT_FIRST_STEP) {
    ws = 0.0;
  }

  cost =
      ws * cost_smooth + wc * cost_colli + wv * cost_vel + wa * cost_acc + 1e-3;

  // cout << "smooth cost:" << ws * cost_smooth << "  collision cost"
  //      << wc * cost_colli << " vel cost " << wv * cost_vel
  //      << "  acc cost: " << wa * cost_acc << "  total:" << cost << endl;

  // sum up all gradient and convert
  gradient = ws * g_smooth + wc * g_colli + wv * g_vel + wa * g_acc;
  _grad.resize(num_dp * 3);

  for (int i = 0; i < num_dp; ++i) {
    _grad[i] = gradient(0, i) + 1e-5;
    _grad[i + num_dp] = gradient(1, i) + 1e-5;
    _grad[i + 2 * num_dp] = gradient(2, i) + 1e-5;
  }

  // get total time
  ros::Time te1 = ros::Time::now();
  total_time += (te1.toSec() - tb1.toSec());

  /* evaluation */
  double time_now = (te1 - time_start).toSec();
  vec_time.push_back(time_now);
  if (vec_cost.size() == 0) {
    vec_cost.push_back(cost);
  } else if (vec_cost.back() > cost) {
    vec_cost.push_back(cost);
  } else {
    vec_cost.push_back(vec_cost.back());
  }
}

// get position from coefficient
void GradTrajOptimizer::getPositionFromCoeff(Eigen::Vector3d &pos,
                                             const Eigen::MatrixXd &coeff,
                                             const int &index,
                                             const double &time) {
  int s = index;
  double t = time;
  float x = coeff(s, 0) + coeff(s, 1) * t + coeff(s, 2) * pow(t, 2) +
            coeff(s, 3) * pow(t, 3) + coeff(s, 4) * pow(t, 4) +
            coeff(s, 5) * pow(t, 5);
  float y = coeff(s, 6) + coeff(s, 7) * t + coeff(s, 8) * pow(t, 2) +
            coeff(s, 9) * pow(t, 3) + coeff(s, 10) * pow(t, 4) +
            coeff(s, 11) * pow(t, 5);
  float z = coeff(s, 12) + coeff(s, 13) * t + coeff(s, 14) * pow(t, 2) +
            coeff(s, 15) * pow(t, 3) + coeff(s, 16) * pow(t, 4) +
            coeff(s, 17) * pow(t, 5);

  pos(0) = x, pos(1) = y, pos(2) = z;
}

// get velocity from cofficient
void GradTrajOptimizer::getVelocityFromCoeff(Eigen::Vector3d &vel,
                                             const Eigen::MatrixXd &coeff,
                                             const int &index,
                                             const double &time) {
  int s = index;
  double t = time;
  float vx = coeff(s, 1) + 2 * coeff(s, 2) * pow(t, 1) +
             3 * coeff(s, 3) * pow(t, 2) + 4 * coeff(s, 4) * pow(t, 3) +
             5 * coeff(s, 5) * pow(t, 4);
  float vy = coeff(s, 7) + 2 * coeff(s, 8) * pow(t, 1) +
             3 * coeff(s, 9) * pow(t, 2) + 4 * coeff(s, 10) * pow(t, 3) +
             5 * coeff(s, 11) * pow(t, 4);
  float vz = coeff(s, 13) + 2 * coeff(s, 14) * pow(t, 1) +
             3 * coeff(s, 15) * pow(t, 2) + 4 * coeff(s, 16) * pow(t, 3) +
             5 * coeff(s, 17) * pow(t, 4);

  vel(0) = vx, vel(1) = vy, vel(2) = vz;
}

// get acceleration from coefficient
void GradTrajOptimizer::getAccelerationFromCoeff(Eigen::Vector3d &acc,
                                                 const Eigen::MatrixXd &coeff,
                                                 const int &index,
                                                 const double &time) {
  int s = index;
  double t = time;
  float ax = 2 * coeff(s, 2) + 6 * coeff(s, 3) * pow(t, 1) +
             12 * coeff(s, 4) * pow(t, 2) + 20 * coeff(s, 5) * pow(t, 3);
  float ay = 2 * coeff(s, 8) + 6 * coeff(s, 9) * pow(t, 1) +
             12 * coeff(s, 10) * pow(t, 2) + 20 * coeff(s, 11) * pow(t, 3);
  float az = 2 * coeff(s, 14) + 6 * coeff(s, 15) * pow(t, 1) +
             12 * coeff(s, 16) * pow(t, 2) + 20 * coeff(s, 17) * pow(t, 3);

  acc(0) = ax, acc(1) = ay, acc(2) = az;
}

inline void GradTrajOptimizer::getDistancePenalty(const double &d,
                                                  double &cost) {
  cost = this->alpha * exp(-(d - this->d0) / this->r);
}

inline void GradTrajOptimizer::getDistancePenaltyGradient(const double &d,
                                                          double &grad) {
  grad = -(this->alpha / this->r) * exp(-(d - this->d0) / this->r);
}

inline void GradTrajOptimizer::getVelocityPenalty(const double &v,
                                                  double &cost) {
  cost = alpha_v * exp((abs(v) - v0) / r_v);
}

inline void GradTrajOptimizer::getVelocityPenaltyGradient(const double &v,
                                                          double &grad) {
  grad = (alpha_v / r_v) * exp((abs(v) - v0) / r_v);
}

inline void GradTrajOptimizer::getAccelerationPenalty(const double &a,
                                                      double &cost) {
  cost = alpha_a * exp((abs(a) - a0) / r_a);
}

inline void GradTrajOptimizer::getAccelerationPenaltyGradient(const double &a,
                                                              double &grad) {
  grad = (alpha_a / r_a) * exp((abs(a) - a0) / r_a);
}

// get distance in signed distance field ,by position query
void GradTrajOptimizer::getDistanceAndGradient(Eigen::Vector3d &pos,
                                               double &dist,
                                               Eigen::Vector3d &grad) {
  dist = sdf_map.getDistWithGradTrilinear(pos, grad);
}

void GradTrajOptimizer::getTimeMatrix(const double &t, Eigen::MatrixXd &T) {
  T.resize(1, 6);
  T.setZero();

  for (int i = 0; i < 6; ++i) {
    T(0, i) = pow(t, i);
  }
}

/** NLopt format cost function */
double GradTrajOptimizer::costFunc(const std::vector<double> &x,
                                   std::vector<double> &grad, void *func_data) {
  GradTrajOptimizer *gtop = reinterpret_cast<GradTrajOptimizer *>(func_data);
  double cost;

  gtop->getCostAndGradient(x, cost, grad);

  return cost;
}