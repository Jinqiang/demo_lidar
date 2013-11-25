#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ros/ros.h>

#include <sensor_msgs/LaserScan.h>

#include "pointDefinition.h"

const double PI = 3.1415926;

const double minAngle = -1.57079637051;
const double increAngle = 0.00436332309619;
const int increNum = 720;

pcl::PointCloud<pcl::PointXYZ>::Ptr syncCloud(new pcl::PointCloud<pcl::PointXYZ>());
ros::Publisher *syncCloudPubPointer = NULL;

void laserDataHandler(const sensor_msgs::LaserScan::ConstPtr& laserData)
{
  pcl::PointXYZ point;
  syncCloud->clear();
  for (int i = 0; i < increNum; i++) {
    double angle = minAngle + increAngle * i;
    point.x = laserData->ranges[i] * sin(angle);
    point.y = 0;
    point.z = laserData->ranges[i] * cos(angle);

    syncCloud->push_back(point);
  }

  syncCloud->header.frame_id = "camera2";
  syncCloud->header.stamp = laserData->header.stamp;
  sensor_msgs::PointCloud2 syncCloud2;
  pcl::toROSMsg(*syncCloud, syncCloud2);
  syncCloudPubPointer->publish(syncCloud2);
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "scanConverter");
  ros::NodeHandle nh;

  ros::Subscriber laserDataSub = nh.subscribe<sensor_msgs::LaserScan> ("/scan", 5, laserDataHandler);

  ros::Publisher syncCloudPub = nh.advertise<sensor_msgs::PointCloud2> ("/sync_scan_cloud_filtered", 5);
  syncCloudPubPointer = &syncCloudPub;

  ros::spin();

  return 0;
}
