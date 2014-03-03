#ifndef point_definition_h
#define point_definition_h

#include <sensor_msgs/PointCloud2.h>
#include <geometry_msgs/PointStamped.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>

struct ImagePoint {
     float u, v;
     int ind;
};

POINT_CLOUD_REGISTER_POINT_STRUCT (ImagePoint,
                                   (float, u, u)
                                   (float, v, v)
                                   (int, ind, ind))

struct DepthPoint {
     float u, v;
     float depth;
     int label;
     int ind;
};

POINT_CLOUD_REGISTER_POINT_STRUCT (DepthPoint,
                                   (float, u, u)
                                   (float, v, v)
                                   (float, depth, depth)
                                   (int, label, label)
                                   (int, ind, ind))

#endif

