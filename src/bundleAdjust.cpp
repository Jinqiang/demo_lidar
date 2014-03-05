#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ros/ros.h>

#include <nav_msgs/Odometry.h>

#include <tf/transform_datatypes.h>
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>

#include "isam/isam.h"
#include "isam/slam_depthmono.h"
#include "isam/robust.h"

#include "cameraParameters.h"
#include "pointDefinition.h"

using namespace std;
using namespace isam;
using namespace Eigen;

const double PI = 3.1415926;

const int keyframeNum = 5;
pcl::PointCloud<DepthPoint>::Ptr depthPoints[keyframeNum];
pcl::PointCloud<DepthPoint>::Ptr depthPointsStacked(new pcl::PointCloud<DepthPoint>());

double depthPointsTime;
bool newKeyframe = false;

double rollRec, pitchRec, yawRec;
double txRec, tyRec, tzRec;

double transformBefBA[6] = {0};
double transformAftBA[6] = {0};

void diffRotation(double cx, double cy, double cz, double lx, double ly, double lz, 
                  double &ox, double &oy, double &oz)
{
  double srx = cos(cx)*cos(cy)*(sin(ly)*sin(lz) + cos(ly)*cos(lz)*sin(lx)) 
             - cos(cx)*sin(cy)*(cos(ly)*sin(lz) - cos(lz)*sin(lx)*sin(ly)) - cos(lx)*cos(lz)*sin(cx);
  ox = -asin(srx);

  double srycrx = cos(cx)*sin(cy)*(cos(ly)*cos(lz) + sin(lx)*sin(ly)*sin(lz)) 
                - cos(cx)*cos(cy)*(cos(lz)*sin(ly) - cos(ly)*sin(lx)*sin(lz)) - cos(lx)*sin(cx)*sin(lz);
  double crycrx = sin(cx)*sin(lx) + cos(cx)*cos(cy)*cos(lx)*cos(ly) + cos(cx)*cos(lx)*sin(cy)*sin(ly);
  oy = atan2(srycrx / cos(ox), crycrx / cos(ox));

  double srzcrx = cos(cx)*cos(lx)*cos(lz)*sin(cz) - (cos(cz)*sin(cy) 
                - cos(cy)*sin(cx)*sin(cz))*(sin(ly)*sin(lz) + cos(ly)*cos(lz)*sin(lx)) 
                - (cos(cy)*cos(cz) + sin(cx)*sin(cy)*sin(cz))*(cos(ly)*sin(lz) - cos(lz)*sin(lx)*sin(ly));
  double crzcrx = (sin(cy)*sin(cz) + cos(cy)*cos(cz)*sin(cx))*(sin(ly)*sin(lz) 
                + cos(ly)*cos(lz)*sin(lx)) + (cos(cy)*sin(cz) - cos(cz)*sin(cx)*sin(cy))*(cos(ly)*sin(lz) 
                - cos(lz)*sin(lx)*sin(ly)) + cos(cx)*cos(cz)*cos(lx)*cos(lz);
  oz = atan2(srzcrx / cos(ox), crzcrx / cos(ox));
}

void transformAssociateToBA()
{
  double txDiff = txRec - transformBefBA[3];
  double tyDiff = tyRec - transformBefBA[4];
  double tzDiff = tzRec - transformBefBA[5];

  double x1 = cos(yawRec) * txDiff - sin(yawRec) * tzDiff;
  double y1 = tyDiff;
  double z1 = sin(yawRec) * txDiff + cos(yawRec) * tzDiff;

  double x2 = x1;
  double y2 = cos(pitchRec) * y1 + sin(pitchRec) * z1;
  double z2 = -sin(pitchRec) * y1 + cos(pitchRec) * z1;

  txDiff = cos(rollRec) * x2 + sin(rollRec) * y2;
  tyDiff = -sin(rollRec) * x2 + cos(rollRec) * y2;
  tzDiff = z2;

  double sbcx = sin(pitchRec);
  double cbcx = cos(pitchRec);
  double sbcy = sin(yawRec);
  double cbcy = cos(yawRec);
  double sbcz = sin(rollRec);
  double cbcz = cos(rollRec);

  double sblx = sin(transformBefBA[0]);
  double cblx = cos(transformBefBA[0]);
  double sbly = sin(transformBefBA[1]);
  double cbly = cos(transformBefBA[1]);
  double sblz = sin(transformBefBA[2]);
  double cblz = cos(transformBefBA[2]);

  double salx = sin(transformAftBA[0]);
  double calx = cos(transformAftBA[0]);
  double saly = sin(transformAftBA[1]);
  double caly = cos(transformAftBA[1]);
  double salz = sin(transformAftBA[2]);
  double calz = cos(transformAftBA[2]);

  double srx = -sbcx*(salx*sblx + calx*caly*cblx*cbly + calx*cblx*saly*sbly) 
             - cbcx*cbcz*(calx*saly*(cbly*sblz - cblz*sblx*sbly) 
             - calx*caly*(sbly*sblz + cbly*cblz*sblx) + cblx*cblz*salx) 
             - cbcx*sbcz*(calx*caly*(cblz*sbly - cbly*sblx*sblz) 
             - calx*saly*(cbly*cblz + sblx*sbly*sblz) + cblx*salx*sblz);
  pitchRec = -asin(srx);

  double srycrx = (cbcy*sbcz - cbcz*sbcx*sbcy)*(calx*saly*(cbly*sblz - cblz*sblx*sbly) 
                - calx*caly*(sbly*sblz + cbly*cblz*sblx) + cblx*cblz*salx) 
                - (cbcy*cbcz + sbcx*sbcy*sbcz)*(calx*caly*(cblz*sbly - cbly*sblx*sblz) 
                - calx*saly*(cbly*cblz + sblx*sbly*sblz) + cblx*salx*sblz) 
                + cbcx*sbcy*(salx*sblx + calx*caly*cblx*cbly + calx*cblx*saly*sbly);
  double crycrx = (cbcz*sbcy - cbcy*sbcx*sbcz)*(calx*caly*(cblz*sbly - cbly*sblx*sblz) 
                - calx*saly*(cbly*cblz + sblx*sbly*sblz) + cblx*salx*sblz) 
                - (sbcy*sbcz + cbcy*cbcz*sbcx)*(calx*saly*(cbly*sblz - cblz*sblx*sbly) 
                - calx*caly*(sbly*sblz + cbly*cblz*sblx) + cblx*cblz*salx) 
                + cbcx*cbcy*(salx*sblx + calx*caly*cblx*cbly + calx*cblx*saly*sbly);
  yawRec = atan2(srycrx / cos(pitchRec), crycrx / cos(pitchRec));
  
  double srzcrx = sbcx*(cblx*cbly*(calz*saly - caly*salx*salz) 
                - cblx*sbly*(caly*calz + salx*saly*salz) + calx*salz*sblx) 
                - cbcx*cbcz*((caly*calz + salx*saly*salz)*(cbly*sblz - cblz*sblx*sbly) 
                + (calz*saly - caly*salx*salz)*(sbly*sblz + cbly*cblz*sblx) 
                - calx*cblx*cblz*salz) + cbcx*sbcz*((caly*calz + salx*saly*salz)*(cbly*cblz 
                + sblx*sbly*sblz) + (calz*saly - caly*salx*salz)*(cblz*sbly - cbly*sblx*sblz) 
                + calx*cblx*salz*sblz);
  double crzcrx = sbcx*(cblx*sbly*(caly*salz - calz*salx*saly) 
                - cblx*cbly*(saly*salz + caly*calz*salx) + calx*calz*sblx) 
                + cbcx*cbcz*((saly*salz + caly*calz*salx)*(sbly*sblz + cbly*cblz*sblx) 
                + (caly*salz - calz*salx*saly)*(cbly*sblz - cblz*sblx*sbly) 
                + calx*calz*cblx*cblz) - cbcx*sbcz*((saly*salz + caly*calz*salx)*(cblz*sbly 
                - cbly*sblx*sblz) + (caly*salz - calz*salx*saly)*(cbly*cblz + sblx*sbly*sblz) 
                - calx*calz*cblx*sblz);
  rollRec = atan2(srzcrx / cos(pitchRec), crzcrx / cos(pitchRec));

  x1 = cos(rollRec) * txDiff - sin(rollRec) * tyDiff;
  y1 = sin(rollRec) * txDiff + cos(rollRec) * tyDiff;
  z1 = tzDiff;

  x2 = x1;
  y2 = cos(pitchRec) * y1 - sin(pitchRec) * z1;
  z2 = sin(pitchRec) * y1 + cos(pitchRec) * z1;

  txDiff = cos(yawRec) * x2 + sin(yawRec) * z2;
  tyDiff = y2;
  tzDiff = -sin(yawRec) * x2 + cos(yawRec) * z2;

  txRec = transformAftBA[3] + txDiff;
  tyRec = transformAftBA[4] + tyDiff;
  tzRec = transformAftBA[5] + tzDiff;
}

void depthPointsHandler(const sensor_msgs::PointCloud2ConstPtr& depthPoints2)
{
  depthPointsTime = depthPoints2->header.stamp.toSec();

  depthPointsStacked->clear();
  pcl::fromROSMsg(*depthPoints2, *depthPointsStacked);
  int depthPointsStackedNum = depthPointsStacked->points.size();

  for (int i = 0; i < keyframeNum; i++) {
    depthPoints[i]->clear();
  }

  int keyframeCount = -1;
  for (int i = 0; i < depthPointsStackedNum; i++) {
    if (depthPointsStacked->points[i].ind == -2) {
      keyframeCount++;
    }

    if (keyframeCount >= 0 && keyframeCount < keyframeNum) {
      depthPoints[keyframeCount]->push_back(depthPointsStacked->points[i]);
    }
  }

  bool enoughPoints = true;
  for (int i = 0; i < keyframeNum; i++) {
    if (depthPoints[i]->points.size() < 10) {
      enoughPoints = false;
    }
  }

  if (enoughPoints) {
    newKeyframe = true;
  }
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "bundleAdjust");
  ros::NodeHandle nh;

  for (int i = 0; i < keyframeNum; i++) {
    pcl::PointCloud<DepthPoint>::Ptr depthPointsTemp(new pcl::PointCloud<DepthPoint>());
    depthPoints[i] = depthPointsTemp;
  }

  ros::Subscriber depthPointsSub = nh.subscribe<sensor_msgs::PointCloud2>
                                   ("/depth_points_stacked", 1, depthPointsHandler);

  ros::Publisher odomBefBAPub = nh.advertise<nav_msgs::Odometry> ("/bef_ba_to_init", 1);

  ros::Publisher odomAftBAPub = nh.advertise<nav_msgs::Odometry> ("/aft_ba_to_init", 1);

  tf::TransformBroadcaster tfBroadcaster;

  Vector2d pp(0, 0);
  DepthmonoCamera camera(1, pp);

  MatrixXd mpNoise = eye(3);
  mpNoise(2, 2) = 0.01;
  Noise noise0 = Information(mpNoise);

  MatrixXd dpNoise = eye(3);
  //dpNoise(2, 2) = 1.;
  Noise noise1 = Information(dpNoise);

  MatrixXd ssNoise = 10000. * eye(6);
  //ssNoise(3, 3) = 100;
  //ssNoise(4, 4) = 100;
  //ssNoise(5, 5) = 100;
  Noise noise2 = Information(ssNoise);

  MatrixXd psNoise = 10000. * eye(6);
  psNoise(3, 3) = 100;
  psNoise(4, 4) = 100;
  psNoise(5, 5) = 100;
  Noise noise3 = Information(psNoise);

  bool status = ros::ok();
  while (status) {
    ros::spinOnce();

    if (newKeyframe) {
      newKeyframe = false;

      Slam ba;

      vector<Pose3d_Node*> poses;
      Pose3d_Node* pose0 = new Pose3d_Node();
      poses.push_back(pose0);
      ba.add_node(pose0);

      rollRec = depthPoints[0]->points[0].depth;
      pitchRec = depthPoints[0]->points[0].u;
      yawRec = depthPoints[0]->points[0].v;
      txRec = depthPoints[0]->points[1].u;
      tyRec = depthPoints[0]->points[1].v;
      tzRec = depthPoints[0]->points[1].depth;

      transformAssociateToBA();

      Pose3d_Factor* poseFactors0 = new Pose3d_Factor(pose0, 
                                    Pose3d(tzRec, txRec, tyRec, yawRec, pitchRec, rollRec), noise2);
      ba.add_factor(poseFactors0);

      rollRec = depthPoints[0]->points[0].depth;
      pitchRec = depthPoints[0]->points[0].u;
      yawRec = depthPoints[0]->points[0].v;
      txRec = depthPoints[0]->points[1].u;
      tyRec = depthPoints[0]->points[1].v;
      tzRec = depthPoints[0]->points[1].depth;

      vector<Pose3d_Pose3d_Factor*> posePoseFactors;
      for (int i = 1; i < keyframeNum; i++) {
        Pose3d_Node* posen = new Pose3d_Node();
        poses.push_back(posen);
        ba.add_node(posen);

        double roll = depthPoints[i]->points[0].depth;
        double pitch = depthPoints[i]->points[0].u;
        double yaw = depthPoints[i]->points[0].v;
        double tx = depthPoints[i]->points[1].u;
        double ty = depthPoints[i]->points[1].v;
        double tz = depthPoints[i]->points[1].depth;

        double txDiff = tx - txRec;
        double tyDiff = ty - tyRec;
        double tzDiff = tz - tzRec;

        double x1 = cos(yawRec) * txDiff - sin(yawRec) * tzDiff;
        double y1 = tyDiff;
        double z1 = sin(yawRec) * txDiff + cos(yawRec) * tzDiff;

        double x2 = x1;
        double y2 = cos(pitchRec) * y1 + sin(pitchRec) * z1;
        double z2 = -sin(pitchRec) * y1 + cos(pitchRec) * z1;

        txDiff = cos(rollRec) * x2 + sin(rollRec) * y2;
        tyDiff = -sin(rollRec) * x2 + cos(rollRec) * y2;
        tzDiff = z2;

        double rollDiff, pitchDiff, yawDiff;
        diffRotation(pitch, yaw, roll, pitchRec, yawRec, rollRec, pitchDiff, yawDiff, rollDiff);

        Pose3d_Pose3d_Factor* poseposeFactorn = new Pose3d_Pose3d_Factor
                                                (poses[i - 1], posen, Pose3d(tzDiff, txDiff, tyDiff, 
                                                yawDiff, pitchDiff, rollDiff), noise3);
        posePoseFactors.push_back(poseposeFactorn);
        ba.add_factor(poseposeFactorn);

        rollRec = roll; pitchRec = pitch; yawRec = yaw;
        txRec = tx; tyRec = ty; tzRec = tz;
      }

      vector<Point3d_Node*> points;
      std::vector<int> pointsInd;
      vector<Depthmono_Factor*> depthmonoFactors;
      for (int i = 0; i < keyframeNum; i++) {
        pcl::PointCloud<DepthPoint>::Ptr dpPointer = depthPoints[i];
        int kfptNum = dpPointer->points.size();
        int ptRecNum = points.size();

        if (i == 0) {
          pcl::PointCloud<DepthPoint>::Ptr dpPointerNext = depthPoints[i + 1];
          int kfptNumNext = dpPointerNext->points.size();

          int ptCountNext = 2;
          for (int j = 2; j < kfptNum; j++) {
            bool ptFound = false;
            for (; ptCountNext < kfptNumNext; ptCountNext++) {
              if (dpPointerNext->points[ptCountNext].ind == dpPointer->points[j].ind) {
                ptFound = true;
              }
              if (dpPointerNext->points[ptCountNext].ind >= dpPointer->points[j].ind) {
                break;
              }
            }

            if (ptFound && dpPointer->points[j].label == 1) {
              Point3d_Node* pointn = new Point3d_Node();
              points.push_back(pointn);
              pointsInd.push_back(dpPointer->points[j].ind);
              ba.add_node(pointn);

              Depthmono_Factor* depthmonoFactorn = new Depthmono_Factor(poses[i], pointn, &camera,
                                                   DepthmonoMeasurement(dpPointer->points[j].u,
                                                   dpPointer->points[j].v, dpPointer->points[j].depth),
                                                   noise1);
              depthmonoFactors.push_back(depthmonoFactorn);
              ba.add_factor(depthmonoFactorn);
            }
          }
        } else if (i == keyframeNum - 1) {
          pcl::PointCloud<DepthPoint>::Ptr dpPointerLast = depthPoints[i - 1];
          int kfptNumLast = dpPointerLast->points.size();

          int ptCountLast = 2;
          int ptRecCount = 0;
          for (int j = 2; j < kfptNum; j++) {
            bool ptFound = false;
            for (; ptCountLast < kfptNumLast; ptCountLast++) {
              if (dpPointerLast->points[ptCountLast].ind == dpPointer->points[j].ind) {
                ptFound = true;
              }
              if (dpPointerLast->points[ptCountLast].ind >= dpPointer->points[j].ind) {
                break;
              }
            }

            if (ptFound /*&& dpPointer->points[j].label == 1*/) {
              bool prFound = false;
              for (; ptRecCount < ptRecNum; ptRecCount++) {
                if (pointsInd[ptRecCount] == dpPointer->points[j].ind) {
                  prFound = true;
                }
                if (pointsInd[ptRecCount] >= dpPointer->points[j].ind) {
                  break;
                }
              }

              Point3d_Node* pointn;
              Depthmono_Factor* depthmonoFactorn;
              if (prFound) {
                pointn = points[ptRecCount];

                depthmonoFactorn = new Depthmono_Factor(poses[i], pointn, &camera,
                                   DepthmonoMeasurement(dpPointer->points[j].u,
                                   dpPointer->points[j].v, dpPointer->points[j].depth), noise0);
                //continue;
              } else {
                pointn = new Point3d_Node();
                points.push_back(pointn);
                pointsInd.push_back(dpPointer->points[j].ind);
                ba.add_node(pointn);

                depthmonoFactorn = new Depthmono_Factor(poses[i], pointn, &camera,
                                   DepthmonoMeasurement(dpPointer->points[j].u,
                                   dpPointer->points[j].v, dpPointer->points[j].depth), noise1);
              }

              depthmonoFactors.push_back(depthmonoFactorn);
              ba.add_factor(depthmonoFactorn);
            }
          }
        } else {
          pcl::PointCloud<DepthPoint>::Ptr dpPointerNext = depthPoints[i + 1];
          pcl::PointCloud<DepthPoint>::Ptr dpPointerLast = depthPoints[i - 1];
          int kfptNumNext = dpPointerNext->points.size();
          int kfptNumLast = dpPointerLast->points.size();

          int ptCountNext = 2;
          int ptCountLast = 2;
          int ptRecCount = 0;
          for (int j = 2; j < kfptNum; j++) {
            bool ptFound = false;
            for (; ptCountNext < kfptNumNext; ptCountNext++) {
              if (dpPointerNext->points[ptCountNext].ind == dpPointer->points[j].ind) {
                ptFound = true;
              }
              if (dpPointerNext->points[ptCountNext].ind >= dpPointer->points[j].ind) {
                break;
              }
            }

            if (!ptFound) {
              for (; ptCountLast < kfptNumLast; ptCountLast++) {
                if (dpPointerLast->points[ptCountLast].ind == dpPointer->points[j].ind) {
                  ptFound = true;
                }
                if (dpPointerLast->points[ptCountLast].ind >= dpPointer->points[j].ind) {
                  break;
                }
              }
            }

            if (ptFound /*&& dpPointer->points[j].label == 1*/) {
              bool prFound = false;
              for (; ptRecCount < ptRecNum; ptRecCount++) {
                if (pointsInd[ptRecCount] == dpPointer->points[j].ind) {
                  prFound = true;
                }
                if (pointsInd[ptRecCount] >= dpPointer->points[j].ind) {
                    break;
                }
              }

              Point3d_Node* pointn;
              Depthmono_Factor* depthmonoFactorn;
              if (prFound) {
                pointn = points[ptRecCount];

                depthmonoFactorn = new Depthmono_Factor(poses[i], pointn, &camera,
                                   DepthmonoMeasurement(dpPointer->points[j].u,
                                   dpPointer->points[j].v, dpPointer->points[j].depth), noise0);
                //continue;
              } else {
                pointn = new Point3d_Node();
                points.push_back(pointn);
                pointsInd.push_back(dpPointer->points[j].ind);
                ba.add_node(pointn);

                depthmonoFactorn = new Depthmono_Factor(poses[i], pointn, &camera,
                                   DepthmonoMeasurement(dpPointer->points[j].u,
                                   dpPointer->points[j].v, dpPointer->points[j].depth), noise1);
              }

              depthmonoFactors.push_back(depthmonoFactorn);
              ba.add_factor(depthmonoFactorn);
            }
          }
        }
      }

      Properties prop = ba.properties();
      prop.method = DOG_LEG;
      ba.set_properties(prop);
      ba.batch_optimization();

      transformBefBA[0] = depthPoints[keyframeNum - 1]->points[0].u;
      transformBefBA[1] = depthPoints[keyframeNum - 1]->points[0].v;
      transformBefBA[2] = depthPoints[keyframeNum - 1]->points[0].depth;
      transformBefBA[3] = depthPoints[keyframeNum - 1]->points[1].u;
      transformBefBA[4] = depthPoints[keyframeNum - 1]->points[1].v;
      transformBefBA[5] = depthPoints[keyframeNum - 1]->points[1].depth;

      transformAftBA[0] = poses[keyframeNum - 1]->value().pitch();
      transformAftBA[1] = poses[keyframeNum - 1]->value().yaw();
      transformAftBA[2] = poses[keyframeNum - 1]->value().roll();
      transformAftBA[3] = poses[keyframeNum - 1]->value().y();
      transformAftBA[4] = poses[keyframeNum - 1]->value().z();
      transformAftBA[5] = poses[keyframeNum - 1]->value().x();

      transformAftBA[0] = (1 - 0.5) * transformAftBA[0] + 0.5 * transformBefBA[0];
      //transformAftBA[1] = (1 - 0.1) * transformAftBA[1] + 0.1 * transformBefBA[1];
      transformAftBA[2] = (1 - 0.5) * transformAftBA[2] + 0.5 * transformBefBA[2];

      int posesNum = poses.size();
      for (int i = 1; i < posesNum; i++) {
        delete poses[i];
      }
      poses.clear();

      delete poseFactors0;

      int posePoseFactorsNum = posePoseFactors.size();
      for (int i = 1; i < posePoseFactorsNum; i++) {
        delete posePoseFactors[i];
      }
      posePoseFactors.clear();

      int pointsNum = points.size();
      for (int i = 1; i < pointsNum; i++) {
        delete points[i];
      }
      points.clear();

      int depthmonoFactorsNum = depthmonoFactors.size();
      for (int i = 1; i < depthmonoFactorsNum; i++) {
        delete depthmonoFactors[i];
      }
      depthmonoFactors.clear();

      geometry_msgs::Quaternion geoQuat = tf::createQuaternionMsgFromRollPitchYaw
                                (transformBefBA[2], -transformBefBA[0], -transformBefBA[1]);

      nav_msgs::Odometry odomBefBA;
      odomBefBA.header.frame_id = "/camera_init";
      odomBefBA.child_frame_id = "/bef_ba";
      odomBefBA.header.stamp = ros::Time().fromSec(depthPointsTime);
      odomBefBA.pose.pose.orientation.x = -geoQuat.y;
      odomBefBA.pose.pose.orientation.y = -geoQuat.z;
      odomBefBA.pose.pose.orientation.z = geoQuat.x;
      odomBefBA.pose.pose.orientation.w = geoQuat.w;
      odomBefBA.pose.pose.position.x = transformBefBA[3];
      odomBefBA.pose.pose.position.y = transformBefBA[4];
      odomBefBA.pose.pose.position.z = transformBefBA[5];
      odomBefBAPub.publish(odomBefBA);

      geoQuat = tf::createQuaternionMsgFromRollPitchYaw
                (transformAftBA[2], -transformAftBA[0], -transformAftBA[1]);

      nav_msgs::Odometry odomAftBA;
      odomAftBA.header.frame_id = "/camera_init";
      odomAftBA.child_frame_id = "/aft_ba";
      odomAftBA.header.stamp = ros::Time().fromSec(depthPointsTime);
      odomAftBA.pose.pose.orientation.x = -geoQuat.y;
      odomAftBA.pose.pose.orientation.y = -geoQuat.z;
      odomAftBA.pose.pose.orientation.z = geoQuat.x;
      odomAftBA.pose.pose.orientation.w = geoQuat.w;
      odomAftBA.pose.pose.position.x = transformAftBA[3];
      odomAftBA.pose.pose.position.y = transformAftBA[4];
      odomAftBA.pose.pose.position.z = transformAftBA[5];
      odomAftBAPub.publish(odomAftBA);

      tf::StampedTransform tfAftBA;
      tfAftBA.frame_id_ = "/camera_init";
      tfAftBA.child_frame_id_ = "/aft_ba";
      tfAftBA.stamp_ = ros::Time().fromSec(depthPointsTime);
      tfAftBA.setRotation(tf::Quaternion(-geoQuat.y, -geoQuat.z, geoQuat.x, geoQuat.w));
      tfAftBA.setOrigin(tf::Vector3(transformAftBA[3], 
                            transformAftBA[4], transformAftBA[5]));
      tfBroadcaster.sendTransform(tfAftBA);
    }

    status = ros::ok();
    cv::waitKey(10);
  }

  return 0;
}
