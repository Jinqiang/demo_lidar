#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ros/ros.h>

#include <nav_msgs/Odometry.h>
#include <sensor_msgs/Imu.h>

#include <tf/transform_datatypes.h>
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>

#include "cameraParameters.h"
#include "pointDefinition.h"

const double PI = 3.1415926;

pcl::PointCloud<ImagePoint>::Ptr imagePointsCur(new pcl::PointCloud<ImagePoint>());
pcl::PointCloud<ImagePoint>::Ptr imagePointsLast(new pcl::PointCloud<ImagePoint>());
pcl::PointCloud<ImagePoint>::Ptr startPointsCur(new pcl::PointCloud<ImagePoint>());
pcl::PointCloud<ImagePoint>::Ptr startPointsLast(new pcl::PointCloud<ImagePoint>());
pcl::PointCloud<pcl::PointXYZHSV>::Ptr startTransCur(new pcl::PointCloud<pcl::PointXYZHSV>());
pcl::PointCloud<pcl::PointXYZHSV>::Ptr startTransLast(new pcl::PointCloud<pcl::PointXYZHSV>());
pcl::PointCloud<pcl::PointXYZHSV>::Ptr ipRelations(new pcl::PointCloud<pcl::PointXYZHSV>());
pcl::PointCloud<pcl::PointXYZHSV>::Ptr ipRelations2(new pcl::PointCloud<pcl::PointXYZHSV>());
pcl::PointCloud<pcl::PointXYZ>::Ptr imagePointsProj(new pcl::PointCloud<pcl::PointXYZ>());
pcl::PointCloud<DepthPoint>::Ptr depthPointsCur(new pcl::PointCloud<DepthPoint>());
pcl::PointCloud<DepthPoint>::Ptr depthPointsLast(new pcl::PointCloud<DepthPoint>());
pcl::PointCloud<DepthPoint>::Ptr depthPointsSend(new pcl::PointCloud<DepthPoint>());

std::vector<int> ipInd;
std::vector<float> ipy2;

double imagePointsCurTime;
double imagePointsLastTime;

int imagePointsCurNum = 0;
int imagePointsLastNum = 0;

int depthPointsCurNum = 0;
int depthPointsLastNum = 0;

pcl::PointCloud<pcl::PointXYZI>::Ptr depthCloud(new pcl::PointCloud<pcl::PointXYZI>());
pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdTree(new pcl::KdTreeFLANN<pcl::PointXYZI>());

double depthCloudTime;
int depthCloudNum = 0;

std::vector<int> pointSearchInd;
std::vector<float> pointSearchSqrDis;

float transformSum[6] = {0};
float angleSum[3] = {0};

pcl::PointXYZ origin(0, 0, 0), xAxis(1, 0, 0), yAxis(0, 1, 0), zAxis(0, 0, 1);

int imuPointerFront = 0;
int imuPointerLast = -1;
const int imuQueLength = 200;
bool imuInited = false;

float imuRollCur = 0, imuPitchCur = 0, imuYawCur = 0;
float imuRollLast = 0, imuPitchLast = 0, imuYawLast = 0;

float imuYawInit = 0;
double imuTime[imuQueLength] = {0};
float imuRoll[imuQueLength] = {0};
float imuPitch[imuQueLength] = {0};
float imuYaw[imuQueLength] = {0};

ros::Publisher *voDataPubPointer = NULL;
tf::TransformBroadcaster *tfBroadcasterPointer = NULL;
ros::Publisher *depthPointsPubPointer = NULL;
ros::Publisher *imagePointsProjPubPointer = NULL;
ros::Publisher *imageShowPubPointer;

const int showDSRate = 3;

IplImage *image;
sensor_msgs::CvBridge bridge;

void transformPoint(pcl::PointXYZ *point, float *r)
{
  double cosrx = cos(r[0]);
  double sinrx = sin(r[0]);
  double cosry = cos(r[1]);
  double sinry = sin(r[1]);
  double cosrz = cos(r[2]);
  double sinrz = sin(r[2]);

  double x1 = cosry * point->x + sinry * point->z;
  double y1 = point->y;
  double z1 = -sinry * point->x + cosry * point->z;

  double x2 = x1;
  double y2 = cosrx * y1 - sinrx * z1;
  double z2 = sinrx * y1 + cosrx * z1;

  point->x = cosrz * x2 - sinrz * y2 + r[3];
  point->y = sinrz * x2 + cosrz * y2 + r[4];
  point->z = z2 + r[5];
}

void transformBack(pcl::PointXYZ *point, float *r)
{
  double cosrx = cos(r[0]);
  double sinrx = sin(r[0]);
  double cosry = cos(r[1]);
  double sinry = sin(r[1]);
  double cosrz = cos(r[2]);
  double sinrz = sin(r[2]);

  double x1 = cosrz * (point->x - r[3]) + sinrz * (point->y - r[4]);
  double y1 = -sinrz * (point->x - r[3]) + cosrz * (point->y - r[4]);
  double z1 = point->z - r[5];

  double x2 = x1;
  double y2 = cosrx * y1 + sinrx * z1;
  double z2 = -sinrx * y1 + cosrx * z1;

  point->x = cosry * x2 - sinry * z2;
  point->y = y2;
  point->z = sinry * x2 + cosry * z2;
}

void imagePointsHandler(const sensor_msgs::PointCloud2ConstPtr& imagePoints2)
{
  imagePointsLastTime = imagePointsCurTime;
  imagePointsCurTime = imagePoints2->header.stamp.toSec();

  imuRollLast = imuRollCur;
  imuPitchLast = imuPitchCur;
  imuYawLast = imuYawCur;

  float transform[6] = {0};
  if (imuPointerLast >= 0) {
    while (imuPointerFront != imuPointerLast) {
      if (imagePointsCurTime < imuTime[imuPointerFront]) {
        break;
      }
      imuPointerFront = (imuPointerFront + 1) % imuQueLength;
    }

    if (imagePointsCurTime > imuTime[imuPointerFront]) {
      imuRollCur = imuRoll[imuPointerFront];
      imuPitchCur = imuPitch[imuPointerFront];
      imuYawCur = imuYaw[imuPointerFront];
    } else {
      int imuPointerBack = (imuPointerFront + imuQueLength - 1) % imuQueLength;
      float ratioFront = (imagePointsCurTime - imuTime[imuPointerBack]) 
                       / (imuTime[imuPointerFront] - imuTime[imuPointerBack]);
      float ratioBack = (imuTime[imuPointerFront] - imagePointsCurTime) 
                      / (imuTime[imuPointerFront] - imuTime[imuPointerBack]);

      imuRollCur = imuRoll[imuPointerFront] * ratioFront + imuRoll[imuPointerBack] * ratioBack;
      imuPitchCur = imuPitch[imuPointerFront] * ratioFront + imuPitch[imuPointerBack] * ratioBack;
      if (imuYaw[imuPointerFront] - imuYaw[imuPointerBack] > PI) {
        imuYawCur = imuYaw[imuPointerFront] * ratioFront + (imuYaw[imuPointerBack] + 2 * PI) * ratioBack;
      } else if (imuYaw[imuPointerFront] - imuYaw[imuPointerBack] < -PI) {
        imuYawCur = imuYaw[imuPointerFront] * ratioFront + (imuYaw[imuPointerBack] - 2 * PI) * ratioBack;
      } else {
        imuYawCur = imuYaw[imuPointerFront] * ratioFront + imuYaw[imuPointerBack] * ratioBack;
      }
    }

    if (imuInited) {
      //transform[0] -= imuPitchCur - imuPitchLast;
      //transform[1] -= imuYawCur - imuYawLast;
      //transform[2] -= imuRollCur - imuRollLast;
    }
  }

  pcl::PointCloud<ImagePoint>::Ptr imagePointsTemp = imagePointsLast;
  imagePointsLast = imagePointsCur;
  imagePointsCur = imagePointsTemp;

  imagePointsCur->clear();
  pcl::fromROSMsg(*imagePoints2, *imagePointsCur);

  imagePointsLastNum = imagePointsCurNum;
  imagePointsCurNum = imagePointsCur->points.size();

  pcl::PointCloud<ImagePoint>::Ptr startPointsTemp = startPointsLast;
  startPointsLast = startPointsCur;
  startPointsCur = startPointsTemp;

  pcl::PointCloud<pcl::PointXYZHSV>::Ptr startTransTemp = startTransLast;
  startTransLast = startTransCur;
  startTransCur = startTransTemp;

  int j = 0;
  pcl::PointXYZI ips;
  pcl::PointXYZHSV ipr;
  ipRelations->clear();
  ipInd.clear();
  for (int i = 0; i < imagePointsLastNum; i++) {
    bool ipFound = false;
    for (; j < imagePointsCurNum; j++) {
      if (imagePointsCur->points[j].ind == imagePointsLast->points[i].ind) {
        ipFound = true;
      }
      if (imagePointsCur->points[j].ind >= imagePointsLast->points[i].ind) {
        break;
      }
    }

    if (ipFound) {
      ipr.x = imagePointsLast->points[i].u;
      ipr.y = imagePointsLast->points[i].v;
      ipr.z = imagePointsCur->points[j].u;
      ipr.h = imagePointsCur->points[j].v;

      ips.x = 10 * ipr.x;
      ips.y = 10 * ipr.y;
      ips.z = 10;
      
      if (depthCloudNum > 10) {
        kdTree->nearestKSearch(ips, 3, pointSearchInd, pointSearchSqrDis);

        double minDepth, maxDepth;
        if (pointSearchSqrDis[0] < 0.5 && pointSearchInd.size() == 3) {
          pcl::PointXYZI depthPoint = depthCloud->points[pointSearchInd[0]];
          double x1 = depthPoint.x * depthPoint.intensity / 10;
          double y1 = depthPoint.y * depthPoint.intensity / 10;
          double z1 = depthPoint.intensity;
          minDepth = z1;
          maxDepth = z1;

          depthPoint = depthCloud->points[pointSearchInd[1]];
          double x2 = depthPoint.x * depthPoint.intensity / 10;
          double y2 = depthPoint.y * depthPoint.intensity / 10;
          double z2 = depthPoint.intensity;
          minDepth = (z2 < minDepth)? z2 : minDepth;
          maxDepth = (z2 > maxDepth)? z2 : maxDepth;

          depthPoint = depthCloud->points[pointSearchInd[2]];
          double x3 = depthPoint.x * depthPoint.intensity / 10;
          double y3 = depthPoint.y * depthPoint.intensity / 10;
          double z3 = depthPoint.intensity;
          minDepth = (z3 < minDepth)? z3 : minDepth;
          maxDepth = (z3 > maxDepth)? z3 : maxDepth;

          double u = ipr.x;
          double v = ipr.y;
          ipr.s = (x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1) 
                / (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2 + u*y1*z2 - u*y2*z1
                - v*x1*z2 + v*x2*z1 - u*y1*z3 + u*y3*z1 + v*x1*z3 - v*x3*z1 + u*y2*z3 
                - u*y3*z2 - v*x2*z3 + v*x3*z2);
          ipr.v = 1;

          if (maxDepth - minDepth > 2) {
            ipr.s = 0;
            ipr.v = 0;
          } else if (ipr.s - maxDepth > 0.2) {
            ipr.s = maxDepth;
          } else if (ipr.s - minDepth < -0.2) {
            ipr.s = minDepth;
          }
        } else {
          ipr.s = 0;
          ipr.v = 0;
        }
      } else {
        ipr.s = 0;
        ipr.v = 0;
      }

      if (fabs(ipr.v) < 0.5) {
        double disX = transformSum[3] - startTransLast->points[i].h;
        double disY = transformSum[4] - startTransLast->points[i].s;
        double disZ = transformSum[5] - startTransLast->points[i].v;

        if (sqrt(disX * disX + disY * disY + disZ * disZ) > 1) {

          float u0 = startPointsLast->points[i].u;
          float v0 = startPointsLast->points[i].v;
          float u1 = ipr.x;
          float v1 = ipr.y;

          float srx0 = sin(-startTransLast->points[i].x);
          float crx0 = cos(-startTransLast->points[i].x);
          float sry0 = sin(-startTransLast->points[i].y);
          float cry0 = cos(-startTransLast->points[i].y);
          float srz0 = sin(-startTransLast->points[i].z);
          float crz0 = cos(-startTransLast->points[i].z);

          float srx1 = sin(-transformSum[0]);
          float crx1 = cos(-transformSum[0]);
          float sry1 = sin(-transformSum[1]);
          float cry1 = cos(-transformSum[1]);
          float srz1 = sin(-transformSum[2]);
          float crz1 = cos(-transformSum[2]);

          float tx0 = -startTransLast->points[i].h;
          float ty0 = -startTransLast->points[i].s;
          float tz0 = -startTransLast->points[i].v;

          float tx1 = -transformSum[3];
          float ty1 = -transformSum[4];
          float tz1 = -transformSum[5];

          float top0 = (srx1*(ty0 - ty1) + crx1*cry1*(tz0 - tz1) - crx1*sry1*(tx0 
                     - tx1))*(v0*((cry0*srz0 + crz0*srx0*sry0)*(cry1*crz1 - srx1*sry1*srz1) + (sry0*srz0 
                     - cry0*crz0*srx0)*(crz1*sry1 + cry1*srx1*srz1) - crx0*crx1*crz0*srz1) 
                     + u0*((crz0*sry0 + cry0*srx0*srz0)*(crz1*sry1 + cry1*srx1*srz1) + (cry0*crz0 
                     - srx0*sry0*srz0)*(cry1*crz1 - srx1*sry1*srz1) + crx0*crx1*srz0*srz1) 
                     + crx0*cry0*(crz1*sry1 + cry1*srx1*srz1) - crx0*sry0*(cry1*crz1 - srx1*sry1*srz1) 
                     - crx1*srx0*srz1) - ((cry1*crz1 - srx1*sry1*srz1)*(tx0 - tx1) + (crz1*sry1 
                     + cry1*srx1*srz1)*(tz0 - tz1) - crx1*srz1*(ty0 - ty1))*(srx0*srx1 
                     + v0*(crx1*cry1*(sry0*srz0 - cry0*crz0*srx0) - crx1*sry1*(cry0*srz0 
                     + crz0*srx0*sry0) + crx0*crz0*srx1) - u0*(crx1*sry1*(cry0*crz0 - srx0*sry0*srz0) 
                     - crx1*cry1*(crz0*sry0 + cry0*srx0*srz0) + crx0*srx1*srz0) + crx0*crx1*sry0*sry1 
                     + crx0*crx1*cry0*cry1);

          float down0 = u1*(srx0*srx1 + v0*(crx1*cry1*(sry0*srz0 - cry0*crz0*srx0) 
                      - crx1*sry1*(cry0*srz0 + crz0*srx0*sry0) + crx0*crz0*srx1) 
                      - u0*(crx1*sry1*(cry0*crz0 - srx0*sry0*srz0) - crx1*cry1*(crz0*sry0 
                      + cry0*srx0*srz0) + crx0*srx1*srz0) + crx0*crx1*sry0*sry1 + crx0*crx1*cry0*cry1) 
                      - v0*((cry0*srz0 + crz0*srx0*sry0)*(cry1*crz1 - srx1*sry1*srz1) + (sry0*srz0 
                      - cry0*crz0*srx0)*(crz1*sry1 + cry1*srx1*srz1) - crx0*crx1*crz0*srz1) 
                      - u0*((crz0*sry0 + cry0*srx0*srz0)*(crz1*sry1 + cry1*srx1*srz1) + (cry0*crz0 
                      - srx0*sry0*srz0)*(cry1*crz1 - srx1*sry1*srz1) + crx0*crx1*srz0*srz1) 
                      - crx0*cry0*(crz1*sry1 + cry1*srx1*srz1) + crx0*sry0*(cry1*crz1 - srx1*sry1*srz1) 
                      + crx1*srx0*srz1;
 
          float top1 = (srx1*(ty0 - ty1) + crx1*cry1*(tz0 - tz1) - crx1*sry1*(tx0 
                     - tx1))*(v0*((cry0*srz0 + crz0*srx0*sry0)*(cry1*srz1 + crz1*srx1*sry1) + (sry0*srz0 
                     - cry0*crz0*srx0)*(sry1*srz1 - cry1*crz1*srx1) + crx0*crx1*crz0*crz1) 
                     + u0*((crz0*sry0 + cry0*srx0*srz0)*(sry1*srz1 - cry1*crz1*srx1) + (cry0*crz0 
                     - srx0*sry0*srz0)*(cry1*srz1 + crz1*srx1*sry1) - crx0*crx1*crz1*srz0) 
                     + crx0*cry0*(sry1*srz1 - cry1*crz1*srx1) - crx0*sry0*(cry1*srz1 + crz1*srx1*sry1) 
                     + crx1*crz1*srx0) - ((cry1*srz1 + crz1*srx1*sry1)*(tx0 - tx1) + (sry1*srz1 
                     - cry1*crz1*srx1)*(tz0 - tz1) + crx1*crz1*(ty0 - ty1))*(srx0*srx1 
                     + v0*(crx1*cry1*(sry0*srz0 - cry0*crz0*srx0) - crx1*sry1*(cry0*srz0 
                     + crz0*srx0*sry0) + crx0*crz0*srx1) - u0*(crx1*sry1*(cry0*crz0 - srx0*sry0*srz0) 
                     - crx1*cry1*(crz0*sry0 + cry0*srx0*srz0) + crx0*srx1*srz0) + crx0*crx1*sry0*sry1 
                     + crx0*crx1*cry0*cry1);

          float down1 = v1*(srx0*srx1 + v0*(crx1*cry1*(sry0*srz0 - cry0*crz0*srx0) 
                      - crx1*sry1*(cry0*srz0 + crz0*srx0*sry0) + crx0*crz0*srx1) 
                      - u0*(crx1*sry1*(cry0*crz0 - srx0*sry0*srz0) - crx1*cry1*(crz0*sry0 
                      + cry0*srx0*srz0) + crx0*srx1*srz0) + crx0*crx1*sry0*sry1 + crx0*crx1*cry0*cry1) 
                      - v0*((cry0*srz0 + crz0*srx0*sry0)*(cry1*srz1 + crz1*srx1*sry1) + (sry0*srz0 
                      - cry0*crz0*srx0)*(sry1*srz1 - cry1*crz1*srx1) + crx0*crx1*crz0*crz1) 
                      - u0*((crz0*sry0 + cry0*srx0*srz0)*(sry1*srz1 - cry1*crz1*srx1) + (cry0*crz0 
                      - srx0*sry0*srz0)*(cry1*srz1 + crz1*srx1*sry1) - crx0*crx1*crz1*srz0) 
                      - crx0*cry0*(sry1*srz1 - cry1*crz1*srx1) + crx0*sry0*(cry1*srz1 + crz1*srx1*sry1) 
                      - crx1*crz1*srx0;

          float depth0 = top0 / down0;
          float depth1 = top1 / down1;
          if (depth0 > 0.5 && depth0 < 100 && depth1 > 0.5 && depth1 < 100) {
            ipr.s = (depth0 + depth1) / 2;
            ipr.v = 2;
          } else if (depth0 > 0.5 && depth0 < 100) {
            ipr.s = depth0;
            ipr.v = 2;
          } else if (depth1 > 0.5 && depth1 < 100) {
            ipr.s = depth1;
            ipr.v = 2;
          }
        }
      }

      ipRelations->push_back(ipr);
      ipInd.push_back(imagePointsLast->points[i].ind);
    }
  }

  int iterNum = 50;
  pcl::PointXYZHSV ipr2, ipr3, ipr4;
  int ipRelationsNum = ipRelations->points.size();
  int ptNumNoDepthRec = 0;
  int ptNumWithDepthRec = 0;
  double meanValueWithDepthRec = 100000;
  for (int iterCount = 0; iterCount < iterNum; iterCount++) {
    ipRelations2->clear();
    ipy2.clear();
    int ptNumNoDepth = 0;
    int ptNumWithDepth = 0;
    double meanValueNoDepth = 0;
    double meanValueWithDepth = 0;
    for (int i = 0; i < ipRelationsNum; i++) {
      ipr = ipRelations->points[i];

      float u0 = ipr.x;
      float v0 = ipr.y;
      float u1 = ipr.z;
      float v1 = ipr.h;

      float srx = sin(transform[0]);
      float crx = cos(transform[0]);
      float sry = sin(transform[1]);
      float cry = cos(transform[1]);
      float srz = sin(transform[2]);
      float crz = cos(transform[2]);
      float tx = transform[3];
      float ty = transform[4];
      float tz = transform[5];

      if (fabs(ipr.v) < 0.5) {

        ipr2.x = v0*(crz*srx*(tx - tz*u1) - crx*(ty*u1 - tx*v1) + srz*srx*(ty - tz*v1)) 
               - u0*(sry*srx*(ty*u1 - tx*v1) + crz*sry*crx*(tx - tz*u1) + sry*srz*crx*(ty - tz*v1)) 
               + cry*srx*(ty*u1 - tx*v1) + cry*crz*crx*(tx - tz*u1) + cry*srz*crx*(ty - tz*v1);

        ipr2.y = u0*((tx - tz*u1)*(srz*sry - crz*srx*cry) - (ty - tz*v1)*(crz*sry + srx*srz*cry) 
               + crx*cry*(ty*u1 - tx*v1)) - (tx - tz*u1)*(srz*cry + crz*srx*sry) 
               + (ty - tz*v1)*(crz*cry - srx*srz*sry) + crx*sry*(ty*u1 - tx*v1);

        ipr2.z = -u0*((tx - tz*u1)*(cry*crz - srx*sry*srz) + (ty - tz*v1)*(cry*srz + srx*sry*crz)) 
               - (tx - tz*u1)*(sry*crz + cry*srx*srz) - (ty - tz*v1)*(sry*srz - cry*srx*crz) 
               - v0*(crx*crz*(ty - tz*v1) - crx*srz*(tx - tz*u1));

        ipr2.h = cry*crz*srx - v0*(crx*crz - srx*v1) - u0*(cry*srz + crz*srx*sry + crx*sry*v1) 
               - sry*srz + crx*cry*v1;

        ipr2.s = crz*sry - v0*(crx*srz + srx*u1) + u0*(cry*crz + crx*sry*u1 - srx*sry*srz) 
               - crx*cry*u1 + cry*srx*srz;

        ipr2.v = u1*(sry*srz - cry*crz*srx) - v1*(crz*sry + cry*srx*srz) + u0*(u1*(cry*srz + crz*srx*sry) 
               - v1*(cry*crz - srx*sry*srz)) + v0*(crx*crz*u1 + crx*srz*v1);

        float y2 = (ty - tz*v1)*(crz*sry + cry*srx*srz) - (tx - tz*u1)*(sry*srz - cry*crz*srx) 
                - v0*(srx*(ty*u1 - tx*v1) + crx*crz*(tx - tz*u1) + crx*srz*(ty - tz*v1)) 
                + u0*((ty - tz*v1)*(cry*crz - srx*sry*srz) - (tx - tz*u1)*(cry*srz + crz*srx*sry) 
                + crx*sry*(ty*u1 - tx*v1)) - crx*cry*(ty*u1 - tx*v1);

        if (ptNumNoDepthRec < 50 || iterCount < 25 || fabs(y2) < 2 * meanValueWithDepthRec / 10000) {
          double scale = 100;
          ipr2.x *= scale;
          ipr2.y *= scale;
          ipr2.z *= scale;
          ipr2.h *= scale;
          ipr2.s *= scale;
          ipr2.v *= scale;
          y2 *= scale;

          ipRelations2->push_back(ipr2);
          ipy2.push_back(y2);

          ptNumNoDepth++;
        } else {
          ipRelations->points[i].v = -1;
        }
      } else if (fabs(ipr.v - 1) < 0.5 || fabs(ipr.v - 2) < 0.5) {

        float d0 = ipr.s;

        ipr3.x = d0*(cry*srz*crx + cry*u1*srx) - d0*u0*(sry*srz*crx + sry*u1*srx) 
               - d0*v0*(u1*crx - srz*srx);

        ipr3.y = d0*(crz*cry + crx*u1*sry - srx*srz*sry) - d0*u0*(crz*sry - crx*u1*cry + srx*srz*cry);

        ipr3.z = -d0*(sry*srz - cry*srx*crz) - d0*u0*(cry*srz + srx*sry*crz) - crx*d0*v0*crz;

        ipr3.h = 1;

        ipr3.s = 0;

        ipr3.v = -u1;

        float y3 = tx - tz*u1 + d0*(crz*sry - crx*cry*u1 + cry*srx*srz) - d0*v0*(crx*srz + srx*u1) 
                 + d0*u0*(cry*crz + crx*sry*u1 - srx*sry*srz);

        ipr4.x = d0*(cry*v1*srx - cry*crz*crx) + d0*u0*(crz*sry*crx - sry*v1*srx) 
               - d0*v0*(crz*srx + v1*crx);

        ipr4.y = d0*(srz*cry + crz*srx*sry + crx*v1*sry) + d0*u0*(crz*srx*cry - srz*sry + crx*v1*cry);

        ipr4.z = d0*(sry*crz + cry*srx*srz) + d0*u0*(cry*crz - srx*sry*srz) - crx*d0*v0*srz;

        ipr4.h = 0;

        ipr4.s = 1;

        ipr4.v = -v1;

        float y4 = ty - tz*v1 - d0*(cry*crz*srx - sry*srz + crx*cry*v1) + d0*v0*(crx*crz - srx*v1) 
                 + d0*u0*(cry*srz + crz*srx*sry + crx*sry*v1);

        if (ptNumWithDepthRec < 50 || iterCount < 25 || 
            sqrt(y3 * y3 + y4 * y4) < 2 * meanValueWithDepthRec) {
          ipRelations2->push_back(ipr3);
          ipy2.push_back(y3);

          ipRelations2->push_back(ipr4);
          ipy2.push_back(y4);

          ptNumWithDepth++;
          meanValueWithDepth += sqrt(y3 * y3 + y4 * y4);
        } else {
          ipRelations->points[i].v = -1;
        }
      }
    }
    meanValueWithDepth /= (ptNumWithDepth + 0.01);
    ptNumNoDepthRec = ptNumNoDepth;
    ptNumWithDepthRec = ptNumWithDepth;
    meanValueWithDepthRec = meanValueWithDepth;

    int ipRelations2Num = ipRelations2->points.size();
    if (ipRelations2Num > 10) {
      cv::Mat matA(ipRelations2Num, 6, CV_32F, cv::Scalar::all(0));
      cv::Mat matAt(6, ipRelations2Num, CV_32F, cv::Scalar::all(0));
      cv::Mat matAtA(6, 6, CV_32F, cv::Scalar::all(0));
      cv::Mat matB(ipRelations2Num, 1, CV_32F, cv::Scalar::all(0));
      cv::Mat matAtB(6, 1, CV_32F, cv::Scalar::all(0));
      cv::Mat matX(6, 1, CV_32F, cv::Scalar::all(0));

      for (int i = 0; i < ipRelations2Num; i++) {
        ipr2 = ipRelations2->points[i];

        matA.at<float>(i, 0) = ipr2.x;
        matA.at<float>(i, 1) = ipr2.y;
        matA.at<float>(i, 2) = ipr2.z;
        matA.at<float>(i, 3) = ipr2.h;
        matA.at<float>(i, 4) = ipr2.s;
        matA.at<float>(i, 5) = ipr2.v;
        matB.at<float>(i, 0) = -0.1 * ipy2[i];
      }
      cv::transpose(matA, matAt);
      matAtA = matAt * matA;
      matAtB = matAt * matB;
      cv::solve(matAtA, matAtB, matX, cv::DECOMP_QR);

      //if (fabs(matX.at<float>(0, 0)) < 0.1 && fabs(matX.at<float>(1, 0)) < 0.1 && 
      //    fabs(matX.at<float>(2, 0)) < 0.1) {
        transform[0] += matX.at<float>(0, 0);
        transform[1] += matX.at<float>(1, 0);
        transform[2] += matX.at<float>(2, 0);
        transform[3] += matX.at<float>(3, 0);
        transform[4] += matX.at<float>(4, 0);
        transform[5] += matX.at<float>(5, 0);
      //}
    }
  }

  if (!imuInited) {
    imuYawInit = imuYawCur;
    transform[0] -= imuPitchCur;
    transform[2] -= imuRollCur;

    imuInited = true;
  }

  transformPoint(&origin, transform);
  transformPoint(&xAxis, transform);
  transformPoint(&yAxis, transform);
  transformPoint(&zAxis, transform);

  /*imagePointsProj->push_back(origin);
  imagePointsProj->push_back(xAxis);
  imagePointsProj->push_back(yAxis);
  imagePointsProj->push_back(zAxis);*/

  float rx = -asin(yAxis.z - origin.z);
  float ry = -atan2(-(xAxis.z - origin.z) / cos(rx), (zAxis.z - origin.z) / cos(rx));
  float rz = -asin(-(yAxis.x - origin.x) / cos(rx));

  if (imuPointerLast >= 0) {
    transformBack(&origin, transform);
    transformBack(&xAxis, transform);
    transformBack(&yAxis, transform);
    transformBack(&zAxis, transform);

    transform[0] -= 0.1 * (imuPitchCur - rx);
    /*if (imuYawCur - imuYawInit - ry > PI) {
      transform[1] -= 0.1 * (imuYawCur - imuYawInit - ry - 2 * PI);
    } else if (imuYawCur - imuYawInit - ry < -PI) {
      transform[1] -= 0.1 * (imuYawCur - imuYawInit - ry + 2 * PI);
    } else {
      transform[1] -= 0.1 * (imuYawCur - imuYawInit - ry);
    }*/
    transform[2] -= 0.1 * (imuRollCur - rz);

    transformPoint(&origin, transform);
    transformPoint(&xAxis, transform);
    transformPoint(&yAxis, transform);
    transformPoint(&zAxis, transform);

    rx = -asin(yAxis.z - origin.z);
    ry = -atan2(-(xAxis.z - origin.z) / cos(rx), (zAxis.z - origin.z) / cos(rx));
    rz = -asin(-(yAxis.x - origin.x) / cos(rx));
  }

  float x1 = cos(rz) * transform[3] - sin(rz) * transform[4];
  float y1 = sin(rz) * transform[3] + cos(rz) * transform[4];
  float z1 = transform[5];

  float x2 = x1;
  float y2 = cos(rx) * y1 - sin(rx) * z1;
  float z2 = sin(rx) * y1 + cos(rx) * z1;

  float tx = transformSum[3] - (cos(ry) * x2 + sin(ry) * z2);
  float ty = transformSum[4] - y2;
  float tz = transformSum[5] - (-sin(ry) * x2 + cos(ry) * z2);

  transformSum[0] = rx;
  transformSum[1] = ry;
  transformSum[2] = rz;
  transformSum[3] = tx;
  transformSum[4] = ty;
  transformSum[5] = tz;

  pcl::PointXYZHSV spc;
  spc.x = transformSum[0];
  spc.y = transformSum[1];
  spc.z = transformSum[2];
  spc.h = transformSum[3];
  spc.s = transformSum[4];
  spc.v = transformSum[5];

  j = 0;
  for (int i = 0; i < imagePointsCurNum; i++) {
    bool ipFound = false;
    for (; j < imagePointsLastNum; j++) {
      if (imagePointsLast->points[j].ind == imagePointsCur->points[i].ind) {
        ipFound = true;
      }
      if (imagePointsLast->points[j].ind >= imagePointsCur->points[i].ind) {
        break;
      }
    }

    if (ipFound) {
      startPointsCur->push_back(startPointsLast->points[j]);
      startTransCur->push_back(startTransLast->points[j]);
    } else {
      startPointsCur->push_back(imagePointsCur->points[i]);
      startTransCur->push_back(spc);
    }
  }
  startPointsLast->clear();
  startTransLast->clear();

  angleSum[0] -= transform[0];
  angleSum[1] -= transform[1];
  angleSum[2] -= transform[2];

  geometry_msgs::Quaternion geoQuat = tf::createQuaternionMsgFromRollPitchYaw(rz, -rx, -ry);

  nav_msgs::Odometry voData;
  voData.header.frame_id = "/camera_init";
  voData.child_frame_id = "/camera";
  voData.header.stamp = imagePoints2->header.stamp;
  voData.pose.pose.orientation.x = -geoQuat.y;
  voData.pose.pose.orientation.y = -geoQuat.z;
  voData.pose.pose.orientation.z = geoQuat.x;
  voData.pose.pose.orientation.w = geoQuat.w;
  voData.pose.pose.position.x = tx;
  voData.pose.pose.position.y = ty;
  voData.pose.pose.position.z = tz;
  voData.twist.twist.angular.x = angleSum[0];
  voData.twist.twist.angular.y = angleSum[1];
  voData.twist.twist.angular.z = angleSum[2];
  voDataPubPointer->publish(voData);

  tf::StampedTransform voTrans;
  voTrans.frame_id_ = "/camera_init";
  voTrans.child_frame_id_ = "/camera";
  voTrans.stamp_ = imagePoints2->header.stamp;
  voTrans.setRotation(tf::Quaternion(-geoQuat.y, -geoQuat.z, geoQuat.x, geoQuat.w));
  voTrans.setOrigin(tf::Vector3(tx, ty, tz));
  tfBroadcasterPointer->sendTransform(voTrans);

  pcl::PointCloud<DepthPoint>::Ptr depthPointsTemp = depthPointsLast;
  depthPointsLast = depthPointsCur;
  depthPointsCur = depthPointsTemp;

  DepthPoint ipd;
  depthPointsCur->clear();

  ipd.u = transformSum[0];
  ipd.v = transformSum[1];
  ipd.depth = transformSum[2];
  ipd.ind = -2;
  depthPointsCur->push_back(ipd);

  ipd.u = transformSum[3];
  ipd.v = transformSum[4];
  ipd.depth = transformSum[5];
  ipd.ind = -1;
  depthPointsCur->push_back(ipd);

  depthPointsLastNum = depthPointsCurNum;
  depthPointsCurNum = 2;

  j = 0;
  pcl::PointXYZ ipp;
  depthPointsSend->clear();
  imagePointsProj->clear();
  for (int i = 0; i < ipRelationsNum; i++) {
    if (fabs(ipRelations->points[i].v - 1) < 0.5 || fabs(ipRelations->points[i].v - 2) < 0.5) {

      ipd.u = ipRelations->points[i].z;
      ipd.v = ipRelations->points[i].h;
      ipd.depth = ipRelations->points[i].s + transform[5];
      ipd.label = ipRelations->points[i].v;
      ipd.ind = ipInd[i];

      depthPointsCur->push_back(ipd);
      depthPointsCurNum++;

      for (; j < depthPointsLastNum; j++) {
        if (depthPointsLast->points[j].ind < ipInd[i]) {
          depthPointsSend->push_back(depthPointsLast->points[j]);
        } else if (depthPointsLast->points[j].ind > ipInd[i]) {
          break;
        }
      }

      ipd.u = ipRelations->points[i].x;
      ipd.v = ipRelations->points[i].y;
      ipd.depth = ipRelations->points[i].s;

      depthPointsSend->push_back(ipd);

      ipp.x = ipRelations->points[i].x * ipRelations->points[i].s;
      ipp.y = ipRelations->points[i].y * ipRelations->points[i].s;
      ipp.z = ipRelations->points[i].s;

      imagePointsProj->push_back(ipp);
    }
  }

  depthPointsSend->header.frame_id = "camera2";
  depthPointsSend->header.stamp = ros::Time().fromSec(imagePointsLastTime);
  sensor_msgs::PointCloud2 depthPoints2;
  pcl::toROSMsg(*depthPointsSend, depthPoints2);
  depthPointsPubPointer->publish(depthPoints2);

  imagePointsProj->header.frame_id = "camera2";
  imagePointsProj->header.stamp = ros::Time().fromSec(imagePointsLastTime);
  sensor_msgs::PointCloud2 imagePointsProj2;
  pcl::toROSMsg(*imagePointsProj, imagePointsProj2);
  imagePointsProjPubPointer->publish(imagePointsProj2);
}

void depthCloudHandler(const sensor_msgs::PointCloud2ConstPtr& depthCloud2)
{
  depthCloudTime = depthCloud2->header.stamp.toSec();

  depthCloud->clear();
  pcl::fromROSMsg(*depthCloud2, *depthCloud);
  depthCloudNum = depthCloud->points.size();

  if (depthCloudNum > 10) {
    for (int i = 0; i < depthCloudNum; i++) {
      depthCloud->points[i].intensity = depthCloud->points[i].z;
      depthCloud->points[i].x *= 10 / depthCloud->points[i].z;
      depthCloud->points[i].y *= 10 / depthCloud->points[i].z;
      depthCloud->points[i].z = 10;
    }

    kdTree->setInputCloud(depthCloud);
  }
}

void imuDataHandler(const sensor_msgs::Imu::ConstPtr& imuData)
{
  double roll, pitch, yaw;
  tf::Quaternion orientation;
  tf::quaternionMsgToTF(imuData->orientation, orientation);
  tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);

  imuPointerLast = (imuPointerLast + 1) % imuQueLength;

  imuTime[imuPointerLast] = imuData->header.stamp.toSec();
  imuRoll[imuPointerLast] = roll;
  imuPitch[imuPointerLast] = pitch;
  imuYaw[imuPointerLast] = yaw;
}

void imageDataHandler(const sensor_msgs::Image::ConstPtr& imageData) 
{
  image = bridge.imgMsgToCv(imageData, "bgr8");

  int ipRelationsNum = ipRelations->points.size();
  for (int i = 0; i < ipRelationsNum; i++) {
    if (fabs(ipRelations->points[i].v) < 0.5) {
      cvCircle(image, cvPoint((kArray[2] - ipRelations->points[i].z * kArray[0]) / showDSRate,
               (kArray[5] - ipRelations->points[i].h * kArray[4]) / showDSRate), 1, CV_RGB(0, 0, 255), 2);
    } else if (fabs(ipRelations->points[i].v - 1) < 0.5) {
      cvCircle(image, cvPoint((kArray[2] - ipRelations->points[i].z * kArray[0]) / showDSRate,
               (kArray[5] - ipRelations->points[i].h * kArray[4]) / showDSRate), 1, CV_RGB(0, 255, 0), 2);
    } else if (fabs(ipRelations->points[i].v - 2) < 0.5) {
      cvCircle(image, cvPoint((kArray[2] - ipRelations->points[i].z * kArray[0]) / showDSRate,
               (kArray[5] - ipRelations->points[i].h * kArray[4]) / showDSRate), 1, CV_RGB(0,255,0), 2);
    } else {
      cvCircle(image, cvPoint((kArray[2] - ipRelations->points[i].z * kArray[0]) / showDSRate,
               (kArray[5] - ipRelations->points[i].h * kArray[4]) / showDSRate), 1, CV_RGB(255, 0, 0), 2);
    }
  }

  sensor_msgs::Image::Ptr imageShow = bridge.cvToImgMsg(image, "bgr8");
  imageShowPubPointer->publish(imageShow);
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "visualOdometry");
  ros::NodeHandle nh;

  ros::Subscriber imagePointsSub = nh.subscribe<sensor_msgs::PointCloud2>
                                   ("/image_points_last", 5, imagePointsHandler);

  ros::Subscriber depthCloudSub = nh.subscribe<sensor_msgs::PointCloud2> 
                                  ("/depth_cloud", 5, depthCloudHandler);

  ros::Subscriber imuDataSub = nh.subscribe<sensor_msgs::Imu> ("/imu/data", 5, imuDataHandler);

  ros::Publisher voDataPub = nh.advertise<nav_msgs::Odometry> ("/cam_to_init", 5);
  voDataPubPointer = &voDataPub;

  tf::TransformBroadcaster tfBroadcaster;
  tfBroadcasterPointer = &tfBroadcaster;

  ros::Publisher depthPointsPub = nh.advertise<sensor_msgs::PointCloud2> ("/depth_points_last", 5);
  depthPointsPubPointer = &depthPointsPub;

  ros::Publisher imagePointsProjPub = nh.advertise<sensor_msgs::PointCloud2> ("/image_points_proj", 1);
  imagePointsProjPubPointer = &imagePointsProjPub;

  ros::Subscriber imageDataSub = nh.subscribe<sensor_msgs::Image>("/image/show", 1, imageDataHandler);

  ros::Publisher imageShowPub = nh.advertise<sensor_msgs::Image>("/image/show_2", 1);
  imageShowPubPointer = &imageShowPub;

  ros::spin();

  return 0;
}
