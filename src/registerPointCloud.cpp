#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ros/ros.h>

#include <nav_msgs/Odometry.h>

#include <tf/transform_datatypes.h>
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>

#include "pointDefinition.h"

const double PI = 3.1415926;

const int keepVoDataNum = 30;
double voDataTime[keepVoDataNum] = {0};
double voRx[keepVoDataNum] = {0};
double voRy[keepVoDataNum] = {0};
double voRz[keepVoDataNum] = {0};
double voTx[keepVoDataNum] = {0};
double voTy[keepVoDataNum] = {0};
double voTz[keepVoDataNum] = {0};
int voDataInd = -1;
int voRegInd = 0;

pcl::PointCloud<pcl::PointXYZ>::Ptr surroundCloud(new pcl::PointCloud<pcl::PointXYZ>());
pcl::PointCloud<pcl::PointXYZ>::Ptr syncCloud(new pcl::PointCloud<pcl::PointXYZ>());
pcl::PointCloud<pcl::PointXYZ>::Ptr tempCloud(new pcl::PointCloud<pcl::PointXYZ>());

int startCount = -1;
const int startSkipNum = 5;

int showCount = -1;
const int showSkipNum = 15;

ros::Publisher *surroundCloudPubPointer = NULL;

void voDataHandler(const nav_msgs::Odometry::ConstPtr& voData)
{
  double time = voData->header.stamp.toSec();

  double rx, ry, rz;
  geometry_msgs::Quaternion geoQuat = voData->pose.pose.orientation;
  tf::Matrix3x3(tf::Quaternion(geoQuat.z, -geoQuat.x, -geoQuat.y, geoQuat.w)).getRPY(rz, rx, ry);

  double tx = voData->pose.pose.position.x;
  double ty = voData->pose.pose.position.y;
  double tz = voData->pose.pose.position.z;

  voDataInd = (voDataInd + 1) % keepVoDataNum;
  voDataTime[voDataInd] = time;
  voRx[voDataInd] = rx;
  voRy[voDataInd] = ry;
  voRz[voDataInd] = rz;
  voTx[voDataInd] = tx;
  voTy[voDataInd] = ty;
  voTz[voDataInd] = tz;
}

void syncCloudHandler(const sensor_msgs::PointCloud2ConstPtr& syncCloud2)
{
  if (startCount < startSkipNum) {
    startCount++;
    return;
  }

  double time = syncCloud2->header.stamp.toSec();

  syncCloud->clear();
  pcl::fromROSMsg(*syncCloud2, *syncCloud);

  double scaleCur = 1;
  double scaleLast = 0;
  int voPreInd = keepVoDataNum - 1;
  if (voDataInd >= 0) {
    while (voDataTime[voRegInd] <= time && voRegInd != voDataInd) {
      voRegInd = (voRegInd + 1) % keepVoDataNum;
    }

    voPreInd = (voRegInd + keepVoDataNum - 1) % keepVoDataNum;
    double voTimePre = voDataTime[voPreInd];
    double voTimeReg = voDataTime[voRegInd];

    if (voTimeReg - voTimePre < 0.5) {
      double scaleLast =  (voTimeReg - time) / (voTimeReg - voTimePre);
      double scaleCur =  (time - voTimePre) / (voTimeReg - voTimePre);
      if (scaleLast > 1) {
        scaleLast = 1;
      } else if (scaleLast < 0) {
        scaleLast = 0;
      }
      if (scaleCur > 1) {
        scaleCur = 1;
      } else if (scaleCur < 0) {
        scaleCur = 0;
      }
    }
  }

  double rx2 = voRx[voRegInd] * scaleCur + voRx[voPreInd] * scaleLast;
  double ry2;
  if (voRy[voRegInd] - voRy[voPreInd] > PI) {
    ry2 = voRy[voRegInd] * scaleCur + (voRy[voPreInd] + 2 * PI) * scaleLast;
  } else if (voRy[voRegInd] - voRy[voPreInd] < -PI) {
    ry2 = voRy[voRegInd] * scaleCur + (voRy[voPreInd] - 2 * PI) * scaleLast;
  } else {
    ry2 = voRy[voRegInd] * scaleCur + voRy[voPreInd] * scaleLast;
  }
  double rz2 = voRz[voRegInd] * scaleCur + voRz[voPreInd] * scaleLast;

  double tx2 = voTx[voRegInd] * scaleCur + voTx[voPreInd] * scaleLast;
  double ty2 = voTy[voRegInd] * scaleCur + voTy[voPreInd] * scaleLast;
  double tz2 = voTz[voRegInd] * scaleCur + voTz[voPreInd] * scaleLast;

  double cosrx2 = cos(rx2);
  double sinrx2 = sin(rx2);
  double cosry2 = cos(ry2);
  double sinry2 = sin(ry2);
  double cosrz2 = cos(rz2);
  double sinrz2 = sin(rz2);

  pcl::PointXYZ point;
  int syncCloudNum = syncCloud->points.size();
  for (int i = 0; i < syncCloudNum; i++) {
    point = syncCloud->points[i];
    double pointDis = sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
    if (pointDis > 0.3 && pointDis < 5) {
      double x1 = cosrz2 * point.x - sinrz2 * point.y;
      double y1 = sinrz2 * point.x + cosrz2 * point.y;
      double z1 = point.z;

      double x2 = x1;
      double y2 = cosrx2 * y1 + sinrx2 * z1;
      double z2 = -sinrx2 * y1 + cosrx2 * z1;

      point.x = cosry2 * x2 - sinry2 * z2 + tx2;
      point.y = y2 + ty2;
      point.z = sinry2 * x2 + cosry2 * z2 + tz2;

      surroundCloud->push_back(point);
    }
  }

  showCount = (showCount + 1) % (showSkipNum + 1);
  if (showCount != showSkipNum) {
    return;
  }

  tempCloud->clear();
  int surroundCloudNum = surroundCloud->points.size();
  for (int i = 0; i < surroundCloudNum; i++) {
    point = surroundCloud->points[i];

    double xDiff = point.x - voTx[voRegInd];
    double yDiff = point.y - voTy[voRegInd];
    double zDiff = point.z - voTz[voRegInd];

    double pointDis = sqrt(xDiff * xDiff + yDiff * yDiff + zDiff * zDiff);
    if (pointDis < 50) {
      tempCloud->push_back(point);
    }
  }

  surroundCloud->clear();
  pcl::VoxelGrid<pcl::PointXYZ> downSizeFilter;
  downSizeFilter.setInputCloud(tempCloud);
  downSizeFilter.setLeafSize(0.2, 0.2, 0.2);
  downSizeFilter.filter(*surroundCloud);

  sensor_msgs::PointCloud2 surroundCloud2;
  pcl::toROSMsg(*surroundCloud, surroundCloud2);
  surroundCloud2.header.frame_id = "/camera_init";
  surroundCloud2.header.stamp = syncCloud2->header.stamp;
  surroundCloudPubPointer->publish(surroundCloud2);
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "registerPointCloud");
  ros::NodeHandle nh;

  ros::Subscriber voDataSub = nh.subscribe<nav_msgs::Odometry> ("/cam2_to_init", 5, voDataHandler);

  ros::Subscriber syncCloudSub = nh.subscribe<sensor_msgs::PointCloud2>
                                 ("/sync_scan_cloud_filtered", 5, syncCloudHandler);

  ros::Publisher surroundCloudPub = nh.advertise<sensor_msgs::PointCloud2> ("/surround_cloud", 1);
  surroundCloudPubPointer = &surroundCloudPub;

  ros::spin();

  return 0;
}
