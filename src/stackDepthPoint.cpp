#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ros/ros.h>

#include "pointDefinition.h"

const double PI = 3.1415926;

const int keyframeNum = 5;
pcl::PointCloud<DepthPoint>::Ptr depthPoints[keyframeNum];
double depthPointsTime[keyframeNum];
int keyframeCount = 0;
int frameCount = 0;

pcl::PointCloud<DepthPoint>::Ptr depthPointsStacked(new pcl::PointCloud<DepthPoint>());
ros::Publisher *depthPointsPubPointer = NULL;

double lastPubTime = 0;

void depthPointsHandler(const sensor_msgs::PointCloud2ConstPtr& depthPoints2)
{
  frameCount = (frameCount + 1) % 5;
  if (frameCount != 0) {
    return;
  }

  pcl::PointCloud<DepthPoint>::Ptr depthPointsCur = depthPoints[0];
  depthPointsCur->clear();
  pcl::fromROSMsg(*depthPoints2, *depthPointsCur);

  for (int i = 0; i < keyframeNum - 1; i++) {
    depthPoints[i] = depthPoints[i + 1];
    depthPointsTime[i] = depthPointsTime[i + 1];
  }
  depthPoints[keyframeNum - 1] = depthPointsCur;
  depthPointsTime[keyframeNum - 1] = depthPoints2->header.stamp.toSec();

  keyframeCount++;
  if (keyframeCount >= keyframeNum && depthPointsTime[0] >= lastPubTime) {
    depthPointsStacked->clear();
    for (int i = 0; i < keyframeNum; i++) {
      *depthPointsStacked += *depthPoints[i];
    }

    depthPointsStacked->header.frame_id = "camera";
    depthPointsStacked->header.stamp = depthPoints2->header.stamp;
    sensor_msgs::PointCloud2 depthPoints2;
    pcl::toROSMsg(*depthPointsStacked, depthPoints2);
    depthPointsPubPointer->publish(depthPoints2);

    lastPubTime = depthPointsTime[keyframeNum - 1];
  }
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "stackDepthPoint");
  ros::NodeHandle nh;

  for (int i = 0; i < keyframeNum; i++) {
    pcl::PointCloud<DepthPoint>::Ptr depthPointsTemp(new pcl::PointCloud<DepthPoint>());
    depthPoints[i] = depthPointsTemp;
  }

  ros::Subscriber depthPointsSub = nh.subscribe<sensor_msgs::PointCloud2>
                                   ("/depth_points_last", 5, depthPointsHandler);

  ros::Publisher depthPointsPub = nh.advertise<sensor_msgs::PointCloud2> ("/depth_points_stacked", 1);
  depthPointsPubPointer = &depthPointsPub;

  ros::spin();

  return 0;
}
