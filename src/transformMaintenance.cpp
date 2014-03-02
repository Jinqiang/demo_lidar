#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ros/ros.h>

#include <nav_msgs/Odometry.h>

#include <tf/transform_datatypes.h>
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>

const double PI = 3.1415926;
const double rad2deg = 180 / PI;
const double deg2rad = PI / 180;

double timeOdomBefBA;
double timeOdomAftBA;

double rollRec, pitchRec, yawRec;
double txRec, tyRec, tzRec;

double transformBefBA[6] = {0};
double transformAftBA[6] = {0};

ros::Publisher *voData2PubPointer = NULL;
tf::TransformBroadcaster *tfBroadcaster2Pointer = NULL;
nav_msgs::Odometry voData2;
tf::StampedTransform voDataTrans2;

void transformAssociateToBA()
{
  double txDiff = txRec - transformBefBA[3];
  double tyDiff = tyRec - transformBefBA[4];
  double tzDiff = tzRec - transformBefBA[5];

  double x1 = cos(yawRec) * txDiff + sin(yawRec) * tzDiff;
  double y1 = tyDiff;
  double z1 = -sin(yawRec) * txDiff + cos(yawRec) * tzDiff;

  double x2 = x1;
  double y2 = cos(pitchRec) * y1 - sin(pitchRec) * z1;
  double z2 = sin(pitchRec) * y1 + cos(pitchRec) * z1;

  txDiff = cos(rollRec) * x2 + sin(rollRec) * y2;
  tyDiff = -sin(rollRec) * x2 + cos(rollRec) * y2;
  tzDiff = z2;

  double sbcx = sin(-pitchRec);
  double cbcx = cos(-pitchRec);
  double sbcy = sin(-yawRec);
  double cbcy = cos(-yawRec);
  double sbcz = sin(rollRec);
  double cbcz = cos(rollRec);

  double sblx = sin(-transformBefBA[0]);
  double cblx = cos(-transformBefBA[0]);
  double sbly = sin(-transformBefBA[1]);
  double cbly = cos(-transformBefBA[1]);
  double sblz = sin(transformBefBA[2]);
  double cblz = cos(transformBefBA[2]);

  double salx = sin(-transformAftBA[0]);
  double calx = cos(-transformAftBA[0]);
  double saly = sin(-transformAftBA[1]);
  double caly = cos(-transformAftBA[1]);
  double salz = sin(transformAftBA[2]);
  double calz = cos(transformAftBA[2]);

  double srx = -sbcx*(salx*sblx + calx*caly*cblx*cbly + calx*cblx*saly*sbly) 
             - cbcx*cbcz*(calx*saly*(cbly*sblz - cblz*sblx*sbly) 
             - calx*caly*(sbly*sblz + cbly*cblz*sblx) + cblx*cblz*salx) 
             - cbcx*sbcz*(calx*caly*(cblz*sbly - cbly*sblx*sblz) 
             - calx*saly*(cbly*cblz + sblx*sbly*sblz) + cblx*salx*sblz);
  pitchRec = asin(srx);

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
  yawRec = -atan2(srycrx / cos(-pitchRec), crycrx / cos(-pitchRec));

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
  rollRec = atan2(srzcrx / cos(-pitchRec), crzcrx / cos(-pitchRec));

  x1 = cos(rollRec) * txDiff - sin(rollRec) * tyDiff;
  y1 = sin(rollRec) * txDiff + cos(rollRec) * tyDiff;
  z1 = tzDiff;

  x2 = x1;
  y2 = cos(pitchRec) * y1 + sin(pitchRec) * z1;
  z2 = -sin(pitchRec) * y1 + cos(pitchRec) * z1;

  txDiff = cos(yawRec) * x2 - sin(yawRec) * z2;
  tyDiff = y2;
  tzDiff = sin(yawRec) * x2 + cos(yawRec) * z2;

  txRec = transformAftBA[3] + txDiff;
  tyRec = transformAftBA[4] + tyDiff;
  tzRec = transformAftBA[5] + tzDiff;
}

void voDataHandler(const nav_msgs::Odometry::ConstPtr& voData)
{
  if (fabs(timeOdomBefBA - timeOdomAftBA) < 0.005) {

    geometry_msgs::Quaternion geoQuat = voData->pose.pose.orientation;
    tf::Matrix3x3(tf::Quaternion(geoQuat.z, -geoQuat.x, -geoQuat.y, geoQuat.w)).getRPY(rollRec, pitchRec, yawRec);

    txRec = voData->pose.pose.position.x;
    tyRec = voData->pose.pose.position.y;
    tzRec = voData->pose.pose.position.z;

    transformAssociateToBA();

    geoQuat = tf::createQuaternionMsgFromRollPitchYaw(rollRec, pitchRec, yawRec);

    voData2.header.stamp = voData->header.stamp;
    voData2.pose.pose.orientation.x = -geoQuat.y;
    voData2.pose.pose.orientation.y = -geoQuat.z;
    voData2.pose.pose.orientation.z = geoQuat.x;
    voData2.pose.pose.orientation.w = geoQuat.w;
    voData2.pose.pose.position.x = txRec;
    voData2.pose.pose.position.y = tyRec;
    voData2.pose.pose.position.z = tzRec;
    voData2PubPointer->publish(voData2);

    voDataTrans2.stamp_ = voData->header.stamp;
    voDataTrans2.setRotation(tf::Quaternion(-geoQuat.y, -geoQuat.z, geoQuat.x, geoQuat.w));
    voDataTrans2.setOrigin(tf::Vector3(txRec, tyRec, tzRec));
    tfBroadcaster2Pointer->sendTransform(voDataTrans2);
  }
}

void odomBefBAHandler(const nav_msgs::Odometry::ConstPtr& odomBefBA)
{
  timeOdomBefBA = odomBefBA->header.stamp.toSec();

  double roll, pitch, yaw;
  geometry_msgs::Quaternion geoQuat = odomBefBA->pose.pose.orientation;
  tf::Matrix3x3(tf::Quaternion(geoQuat.z, -geoQuat.x, -geoQuat.y, geoQuat.w)).getRPY(roll, pitch, yaw);

  transformBefBA[0] = pitch;
  transformBefBA[1] = yaw;
  transformBefBA[2] = roll;

  transformBefBA[3] = odomBefBA->pose.pose.position.x;
  transformBefBA[4] = odomBefBA->pose.pose.position.y;
  transformBefBA[5] = odomBefBA->pose.pose.position.z;
}

void odomAftBAHandler(const nav_msgs::Odometry::ConstPtr& odomAftBA)
{
  timeOdomAftBA = odomAftBA->header.stamp.toSec();

  double roll, pitch, yaw;
  geometry_msgs::Quaternion geoQuat = odomAftBA->pose.pose.orientation;
  tf::Matrix3x3(tf::Quaternion(geoQuat.z, -geoQuat.x, -geoQuat.y, geoQuat.w)).getRPY(roll, pitch, yaw);

  transformAftBA[0] = pitch;
  transformAftBA[1] = yaw;
  transformAftBA[2] = roll;

  transformAftBA[3] = odomAftBA->pose.pose.position.x;
  transformAftBA[4] = odomAftBA->pose.pose.position.y;
  transformAftBA[5] = odomAftBA->pose.pose.position.z;
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "transformMaintenance");
  ros::NodeHandle nh;

  ros::Subscriber voDataSub = nh.subscribe<nav_msgs::Odometry> 
                              ("/cam_to_init", 1, voDataHandler);

  ros::Subscriber odomBefBASub = nh.subscribe<nav_msgs::Odometry> 
                                 ("/bef_ba_to_init", 1, odomBefBAHandler);

  ros::Subscriber odomAftBASub = nh.subscribe<nav_msgs::Odometry> 
                                 ("/aft_ba_to_init", 1, odomAftBAHandler);

  ros::Publisher voData2Pub = nh.advertise<nav_msgs::Odometry> ("/cam2_to_init", 1);
  voData2PubPointer = &voData2Pub;
  voData2.header.frame_id = "/camera_init";
  voData2.child_frame_id = "/camera2";

  tf::TransformBroadcaster tfBroadcaster2;
  tfBroadcaster2Pointer = &tfBroadcaster2;
  voDataTrans2.frame_id_ = "/camera_init";
  voDataTrans2.child_frame_id_ = "/camera2";

  ros::spin();

  return 0;
}
