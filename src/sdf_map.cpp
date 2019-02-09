#include "grad_traj_optimization/sdf_map.h"

SDFMap::SDFMap(Eigen::Vector3d origin, double resolution,
               Eigen::Vector3d map_size) {
  this->origin = origin;
  this->resolution = resolution;
  this->resolution_inv = 1 / resolution;
  this->map_size = map_size;
  for (int i = 0; i < 3; ++i) grid_size(i) = ceil(map_size(i) / resolution);
  // SETY << "grid num:" << grid_size.transpose() << REC;
  min_range = origin;
  max_range = origin + map_size;
  min_vec = Eigen::Vector3i::Zero();
  max_vec = grid_size - Eigen::Vector3i::Ones();

  // initialize size of buffer
  occupancy_buffer.resize(grid_size(0) * grid_size(1) * grid_size(2));
  distance_buffer.resize(grid_size(0) * grid_size(1) * grid_size(2));
  tmp_buffer1.resize(grid_size(0) * grid_size(1) * grid_size(2));
  tmp_buffer2.resize(grid_size(0) * grid_size(1) * grid_size(2));

  fill(distance_buffer.begin(), distance_buffer.end(), 10000);
  fill(occupancy_buffer.begin(), occupancy_buffer.end(), 0.0);
}

void SDFMap::resetBuffer() { this->resetBuffer(min_range, max_range); }

void SDFMap::resetBuffer(Eigen::Vector3d min_pos, Eigen::Vector3d max_pos) {
  min_pos(0) = max(min_pos(0), min_range(0));
  min_pos(1) = max(min_pos(1), min_range(1));
  min_pos(2) = max(min_pos(2), min_range(2));

  max_pos(0) = min(max_pos(0), max_range(0));
  max_pos(1) = min(max_pos(1), max_range(1));
  max_pos(2) = min(max_pos(2), max_range(2));

  Eigen::Vector3i min_id, max_id;

  posToIndex(min_pos, min_id);
  posToIndex(
      max_pos - Eigen::Vector3d(resolution / 2, resolution / 2, resolution / 2),
      max_id);

  /* reset occ and dist buffer */
  for (int x = min_id(0); x <= max_id(0); ++x)
    for (int y = min_id(1); y <= max_id(1); ++y)
      for (int z = min_id(2); z <= max_id(2); ++z) {
        occupancy_buffer[x * grid_size(1) * grid_size(2) + y * grid_size(2) +
                         z] = 0.0;
        distance_buffer[x * grid_size(1) * grid_size(2) + y * grid_size(2) +
                        z] = 10000;
      }
}

bool SDFMap::isInMap(Eigen::Vector3d pos) {
  if (pos(0) < min_range(0) + 1e-4 || pos(1) < min_range(1) + 1e-4 ||
      pos(2) < min_range(2) + 1e-4) {
    // SETY << "less than min range!" << REC;
    return false;
  }

  if (pos(0) > max_range(0) - 1e-4 || pos(1) > max_range(1) - 1e-4 ||
      pos(2) > max_range(2) - 1e-4) {
    // SETY << "larger than max range!" << REC;
    return false;
  }

  return true;
}

void SDFMap::posToIndex(Eigen::Vector3d pos, Eigen::Vector3i& id) {
  for (int i = 0; i < 3; ++i)
    id(i) = floor((pos(i) - origin(i)) * resolution_inv);
}

void SDFMap::indexToPos(Eigen::Vector3i id, Eigen::Vector3d& pos) {
  for (int i = 0; i < 3; ++i) pos(i) = (id(i) + 0.5) * resolution + origin(i);
}

void SDFMap::setOccupancy(Eigen::Vector3d pos, int occ) {
  if (occ != 1 && occ != 0) {
    // SETY << "occ value error!" << REC;
    cout << "occ value error!" << endl;
    return;
  }

  if (!isInMap(pos)) return;

  Eigen::Vector3i id;
  posToIndex(pos, id);

  // (x, y, z) -> x*ny*nz + y*nz + z
  // SETY << "..."
  //      << id(0) * grid_size(1) * grid_size(2) + id(1) * grid_size(2) + id(2)
  //      << REC;
  // SETY << "..." << occupancy_buffer.size() << REC;
  occupancy_buffer[id(0) * grid_size(1) * grid_size(2) + id(1) * grid_size(2) +
                   id(2)] = occ;
}

int SDFMap::getOccupancy(Eigen::Vector3d pos) {
  if (!isInMap(pos)) return -1;

  Eigen::Vector3i id;
  posToIndex(pos, id);

  // (x, y, z) -> x*ny*nz + y*nz + z
  return occupancy_buffer[id(0) * grid_size(1) * grid_size(2) +
                          id(1) * grid_size(2) + id(2)];
}

int SDFMap::getOccupancy(Eigen::Vector3i id) {
  if (id(0) < 0 || id(0) >= grid_size(0) || id(1) < 0 ||
      id(1) >= grid_size(1) || id(2) < 0 || id(2) >= grid_size(2))
    return -1;

  // (x, y, z) -> x*ny*nz + y*nz + z
  return occupancy_buffer[id(0) * grid_size(1) * grid_size(2) +
                          id(1) * grid_size(2) + id(2)];
}

void SDFMap::getOccupancyMarker(visualization_msgs::Marker& m, int id,
                                Eigen::Vector4d color) {
  m.header.frame_id = "world";
  m.id = id;
  m.type = visualization_msgs::Marker::CUBE_LIST;
  m.action = visualization_msgs::Marker::MODIFY;
  m.scale.x = resolution * 0.9;
  m.scale.y = resolution * 0.9;
  m.scale.z = resolution * 0.9;
  m.color.a = color(3);
  m.color.r = color(0);
  m.color.g = color(1);
  m.color.b = color(2);

  // iterate the map
  for (int x = 0; x < grid_size(0); ++x)
    for (int y = 0; y < grid_size(1); ++y)
      for (int z = 0; z < grid_size(2); ++z) {
        if (1 != occupancy_buffer[x * grid_size(1) * grid_size(2) +
                                  y * grid_size(2) + z])
          continue;

        Eigen::Vector3d pos;
        indexToPos(Eigen::Vector3i(x, y, z), pos);

        geometry_msgs::Point p;
        p.x = pos(0);
        p.y = pos(1);
        p.z = pos(2);
        m.points.push_back(p);
      }
}

double SDFMap::getDistance(Eigen::Vector3d pos) {
  if (!isInMap(pos)) return -1;

  Eigen::Vector3i id;
  posToIndex(pos, id);

  // (x, y, z) -> x*ny*nz + y*nz + z
  return distance_buffer[id(0) * grid_size(1) * grid_size(2) +
                         id(1) * grid_size(2) + id(2)];
}

double SDFMap::getDistance(Eigen::Vector3i id) {
  id(0) = max(min(id(0), grid_size(0) - 1), 0);
  id(1) = max(min(id(1), grid_size(1) - 1), 0);
  id(2) = max(min(id(2), grid_size(2) - 1), 0);

  // (x, y, z) -> x*ny*nz + y*nz + z
  return distance_buffer[id(0) * grid_size(1) * grid_size(2) +
                         id(1) * grid_size(2) + id(2)];
}

double SDFMap::getDistance(int x, int y, int z) {
  x = max(min(x, grid_size(0) - 1), 0);
  y = max(min(y, grid_size(1) - 1), 0);
  z = max(min(z, grid_size(2) - 1), 0);

  return distance_buffer[x * grid_size(1) * grid_size(2) + y * grid_size(2) +
                         z];
}

double SDFMap::getDistWithGradTrilinear(Eigen::Vector3d pos,
                                        Eigen::Vector3d& grad) {
  if (!isInMap(pos)) return -1;

  /* this is too young too simple */
  // gradient(0) = (getDistance(id(0) + 1, id(1), id(2)) -
  //                getDistance(id(0) - 1, id(1), id(2))) *
  //               0.5 * resolution_inv;
  // gradient(1) = (getDistance(id(0), id(1) + 1, id(2)) -
  //                getDistance(id(0), id(1) - 1, id(2))) *
  //               0.5 * resolution_inv;
  // gradient(2) = (getDistance(id(0), id(1), id(2) + 1) -
  //                getDistance(id(0), id(1), id(2) - 1)) *
  //               0.5 * resolution_inv;

  /* use trilinear interpolation */
  Eigen::Vector3d pos_m = pos - 0.5 * resolution * Eigen::Vector3d::Ones();

  Eigen::Vector3i idx;
  posToIndex(pos_m, idx);

  Eigen::Vector3d idx_pos, diff;
  indexToPos(idx, idx_pos);

  diff = (pos - idx_pos) * resolution_inv;

  double values[2][2][2];
  for (int x = 0; x < 2; x++) {
    for (int y = 0; y < 2; y++) {
      for (int z = 0; z < 2; z++) {
        Eigen::Vector3i current_idx = idx + Eigen::Vector3i(x, y, z);
        values[x][y][z] = getDistance(current_idx);
      }
    }
  }

  double v00 = (1 - diff[0]) * values[0][0][0] + diff[0] * values[1][0][0];
  double v01 = (1 - diff[0]) * values[0][0][1] + diff[0] * values[1][0][1];
  double v10 = (1 - diff[0]) * values[0][1][0] + diff[0] * values[1][1][0];
  double v11 = (1 - diff[0]) * values[0][1][1] + diff[0] * values[1][1][1];

  double v0 = (1 - diff[1]) * v00 + diff[1] * v10;
  double v1 = (1 - diff[1]) * v01 + diff[1] * v11;

  double dist = (1 - diff[2]) * v0 + diff[2] * v1;

  grad[2] = (v1 - v0) * resolution_inv;
  grad[1] =
      ((1 - diff[2]) * (v10 - v00) + diff[2] * (v11 - v01)) * resolution_inv;
  grad[0] = (1 - diff[2]) * (1 - diff[1]) * (values[1][0][0] - values[0][0][0]);
  grad[0] += (1 - diff[2]) * diff[1] * (values[1][1][0] - values[0][1][0]);
  grad[0] += diff[2] * (1 - diff[1]) * (values[1][0][1] - values[0][0][1]);
  grad[0] += diff[2] * diff[1] * (values[1][1][1] - values[0][1][1]);

  grad[0] *= resolution_inv;

  return dist;
}

void SDFMap::setUpdateRange(Eigen::Vector3d min_pos, Eigen::Vector3d max_pos) {
  /* chou gou shi! */
  // if (!isInMap(min_pos)) min_pos = min_range;
  // if (!isInMap(max_pos)) max_pos = max_range;
  min_pos(0) = max(min_pos(0), min_range(0));
  min_pos(1) = max(min_pos(1), min_range(1));
  min_pos(2) = max(min_pos(2), min_range(2));

  max_pos(0) = min(max_pos(0), max_range(0));
  max_pos(1) = min(max_pos(1), max_range(1));
  max_pos(2) = min(max_pos(2), max_range(2));

  posToIndex(min_pos, min_vec);
  posToIndex(
      max_pos - Eigen::Vector3d(resolution / 2, resolution / 2, resolution / 2),
      max_vec);
  // SETY << "mind:" << min_pos.transpose() << ", maxd:" << max_pos.transpose()
  //      << REC;
  // SETY << "min:" << min_vec.transpose() << ", max:" << max_vec.transpose()
  //      << REC;
}

template <typename F_get_val, typename F_set_val>
void SDFMap::fillESDF(F_get_val f_get_val, F_set_val f_set_val, int start,
                      int end, int dim) {
  int v[grid_size(dim)];
  double z[grid_size(dim) + 1];

  int k = start;
  v[start] = start;
  z[start] = -std::numeric_limits<double>::max();
  z[start + 1] = std::numeric_limits<double>::max();

  for (int q = start + 1; q <= end; q++) {
    k++;
    double s;

    do {
      k--;
      s = ((f_get_val(q) + q * q) - (f_get_val(v[k]) + v[k] * v[k])) /
          (2 * q - 2 * v[k]);
      // ROS_INFO_STREAM("k: " << k << " s: " <<  s  << " z[k] " << z[k] << "
      // v[k] " << v[k]);

    } while (s <= z[k]);

    k++;

    v[k] = q;
    z[k] = s;
    z[k + 1] = std::numeric_limits<double>::max();
  }

  k = start;

  for (int q = start; q <= end; q++) {
    while (z[k + 1] < q) k++;
    double val = (q - v[k]) * (q - v[k]) + f_get_val(v[k]);
    //      if(val < std::numeric_limits<_Scalar>::max())
    //  ROS_INFO_STREAM("val: " << val << " q: " << q << " v[k] " << v[k]);
    // if(val > truncation_distance_*truncation_distance_) val =
    // std::numeric_limits<_Scalar>::max();
    f_set_val(q, val);
  }
}

void SDFMap::updateESDF3d() {
  for (int x = min_vec[0]; x <= max_vec[0]; x++) {
    for (int y = min_vec[1]; y <= max_vec[1]; y++) {
      fillESDF(
          [&](int z) {
            return occupancy_buffer[x * grid_size(1) * grid_size(2) +
                                    y * grid_size(2) + z] == 1
                       ? 0
                       : std::numeric_limits<double>::max();
          },
          [&](int z, double val) {
            tmp_buffer1[x * grid_size(1) * grid_size(2) + y * grid_size(2) +
                        z] = val;
          },
          min_vec[2], max_vec[2], 2);
    }
  }

  for (int x = min_vec[0]; x <= max_vec[0]; x++) {
    for (int z = min_vec[2]; z <= max_vec[2]; z++) {
      fillESDF(
          [&](int y) {
            // SETY << "get xyz:" << x << ", " << y << ", " << z << REC;
            return tmp_buffer1[x * grid_size(1) * grid_size(2) +
                               y * grid_size(2) + z];
          },
          [&](int y, double val) {
            // SETY << "set xyz:" << x << ", " << y << ", " << z << REC;
            // SETY << "index:" << x * grid_size(1) * grid_size(2) + y *
            // grid_size(2) + z << REC; SETY << "buffer length:" <<
            // tmp_buffer2.size() << REC;
            tmp_buffer2[x * grid_size(1) * grid_size(2) + y * grid_size(2) +
                        z] = val;
          },
          min_vec[1], max_vec[1], 1);
    }
  }

  for (int y = min_vec[1]; y <= max_vec[1]; y++) {
    for (int z = min_vec[2]; z <= max_vec[2]; z++) {
      fillESDF(
          [&](int x) {
            return tmp_buffer2[x * grid_size(1) * grid_size(2) +
                               y * grid_size(2) + z];
          },
          [&](int x, double val) {
            distance_buffer[x * grid_size(1) * grid_size(2) + y * grid_size(2) +
                            z] =
                min(resolution * std::sqrt(val),
                    distance_buffer[x * grid_size(1) * grid_size(2) +
                                    y * grid_size(2) + z]);
          },
          min_vec[0], max_vec[0], 0);
    }
  }

  min_vec = Eigen::Vector3i::Zero();
  max_vec = grid_size - Eigen::Vector3i::Ones();
}

void SDFMap::getESDFMarker(vector<visualization_msgs::Marker>& markers, int id,
                           Eigen::Vector3d color) {
  double max_dist = getMaxDistance();

  // get marker in several distance level
  const int level = ceil(max_dist * resolution_inv);

  for (int i = 0; i < level; ++i) {
    visualization_msgs::Marker m;
    m.header.frame_id = "world";
    m.id = i + level * id;
    m.type = visualization_msgs::Marker::CUBE_LIST;
    m.action = visualization_msgs::Marker::ADD;
    m.scale.x = resolution * 0.9;
    m.scale.y = resolution * 0.9;
    m.scale.z = resolution * 0.9;
    m.color.r = color(0);
    m.color.g = color(1);
    m.color.b = color(2);

    // transparency and distance conversion
    double min_a = 0.05, max_a = 0.25;
    double da = (max_a - min_a) / (level - 1);
    m.color.a = max_a - da * i;
    // SETY << "alpha:" << m.color.a << REC;

    // distance level
    double delta_d = max_dist / level;
    double min_d = i * delta_d - 1e-3;
    double max_d = (i + 1) * delta_d - 1e-3;

    // iterate the map
    for (int x = 0; x < grid_size(0); ++x)
      for (int y = 0; y < grid_size(1); ++y)
        for (int z = 0; z < grid_size(2) - 15; ++z) {
          double dist = distance_buffer[x * grid_size(1) * grid_size(2) +
                                        y * grid_size(2) + z];
          bool in_range = dist < max_d && dist >= min_d;
          if (!in_range) continue;

          Eigen::Vector3d pos;
          indexToPos(Eigen::Vector3i(x, y, z), pos);

          geometry_msgs::Point p;
          p.x = pos(0);
          p.y = pos(1);
          p.z = pos(2);
          m.points.push_back(p);
        }
    markers.push_back(m);
  }
}

double SDFMap::getMaxDistance() {
  // get the max distance
  double max_dist = -1;
  for (int i = 0; i < int(distance_buffer.size()); ++i) {
    if (distance_buffer[i] > max_dist) max_dist = distance_buffer[i];
  }
  // SETY << "Max distance is:" << max_dist << REC;
  return max_dist;
}