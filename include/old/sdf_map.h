#ifndef _SDF_MAP_H
#define _SDF_MAP_H

#include <visualization_msgs/Marker.h>
#include <Eigen/Eigen>
#include <iostream>

using namespace std;

class SDFMap {
 private:
  // data are saved in vector
  std::vector<int> occupancy_buffer;  // 0 is free, 1 is occupied
  std::vector<double> distance_buffer;
  std::vector<double> tmp_buffer1, tmp_buffer2;

  // map property
  Eigen::Vector3d origin, map_size;
  Eigen::Vector3d min_range, max_range;  // map range in pos
  Eigen::Vector3i grid_size;             // map range in index
  double resolution, resolution_inv;
  Eigen::Vector3i min_vec, max_vec;  // the min and max updated range, unit is 1
  double truncated_distance = 20.0;

  bool isInMap(Eigen::Vector3d pos);
  void posToIndex(Eigen::Vector3d pos, Eigen::Vector3i& id);
  void indexToPos(Eigen::Vector3i id, Eigen::Vector3d& pos);

  template <typename F_get_val, typename F_set_val>
  void fillESDF(F_get_val f_get_val, F_set_val f_set_val, int start, int end,
                int dim);

 public:
  SDFMap() {}
  SDFMap(Eigen::Vector3d origin, double resolution, Eigen::Vector3d map_size);
  ~SDFMap() {}

  // occupancy management
  void resetBuffer();
  void resetBuffer(Eigen::Vector3d min, Eigen::Vector3d max);
  void setOccupancy(Eigen::Vector3d pos, int occ = 1);
  int getOccupancy(Eigen::Vector3d pos);
  int getOccupancy(Eigen::Vector3i id);
  void getOccupancyMarker(visualization_msgs::Marker& m, int id,
                          Eigen::Vector4d color);

  // distance field management
  double getDistance(Eigen::Vector3d pos);
  double getDistance(Eigen::Vector3i id);
  double getDistance(int x, int y, int z);
  double getDistWithGradTrilinear(Eigen::Vector3d pos, Eigen::Vector3d& grad);
  void setUpdateRange(Eigen::Vector3d min_pos, Eigen::Vector3d max_pos);
  void updateESDF3d();
  void getESDFMarker(vector<visualization_msgs::Marker>& markers, int id,
                     Eigen::Vector3d color);
  double getMaxDistance();
};

#endif