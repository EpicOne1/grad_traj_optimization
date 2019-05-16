#ifndef _EDT_ENVIRONMENT_H_
#define _EDT_ENVIRONMENT_H_

#include <Eigen/Eigen>
#include <iostream>
#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>

#include <grad_traj_optimization/sdf_map.h>
#include <grad_traj_optimization/obj_predictor.h>

using std::cout;
using std::endl;
using std::list;
using std::shared_ptr;
using std::unique_ptr;
using std::vector;

namespace dyn_planner
{
class EDTEnvironment
{
private:
  /* data */
  SDFMap::Ptr sdf_map_;
  ObjPrediction obj_prediction_;
  ObjScale obj_scale_;
  double resolution_inv_;

  double distToBox(int idx, const Eigen::Vector3d& pos, const double& time);
  double minDistToAllBox(const Eigen::Vector3d& pos, const double& time);

public:
  EDTEnvironment(/* args */) {}
  ~EDTEnvironment() {}

  void init();
  void setMap(SDFMap::Ptr map);
  void setObjPrediction(ObjPrediction prediction);
  void setObjScale(ObjScale scale);

  void evaluateEDTWithGrad(const Eigen::Vector3d& pos, const double& time, double& dist,
                           Eigen::Vector3d& grad);

  double evaluateCoarseEDT(const Eigen::Vector3d& pos, const double& time);

  bool odomValid() { return sdf_map_->odomValid(); }
  bool mapValid() { return sdf_map_->mapValid(); }
  nav_msgs::Odometry getOdom() { return sdf_map_->getOdom(); }
  void getMapRegion(Eigen::Vector3d& ori, Eigen::Vector3d& size) { sdf_map_->getRegion(ori, size); }

  typedef shared_ptr<EDTEnvironment> Ptr;
};

}  // namespace dyn_planner

#endif