#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ros/ros.h>

#include "cameraParameters.h"
#include "opencv2/ocl/ocl.hpp"
#include "pointDefinition.h"

using namespace std;
using namespace cv;
using namespace cv::ocl;

bool systemInited = false;
bool isOddFrame = true;

double timeCur, timeLast;

int detectionCount = 0;
const int detectionSkipNum = 3;

int showCount = 0;
const int showSkipNum = 2;
const int showDSRate = 2;
Size showSize = Size(imageWidth / showDSRate, imageHeight / showDSRate);

Mat image0, image1;
oclMat oclImage0, oclImage1;
Mat imageShow, harrisLast;
oclMat oclImageShow, oclHarrisLast;

Mat kMat = Mat(3, 3, CV_64FC1, kImage);
Mat dMat = Mat(4, 1, CV_64FC1, dImage);

Mat mapx, mapy;

const int maxFeatureNumPerSubregion = 20;
const int xSubregionNum = 6;
const int ySubregionNum = 4;
const int totalSubregionNum = xSubregionNum * ySubregionNum;
const int MAXFEATURENUM = maxFeatureNumPerSubregion * totalSubregionNum;

const int xBoundary = 20;
const int yBoundary = 20;
const double subregionWidth = (double)(imageWidth - 2 * xBoundary) / (double)xSubregionNum;
const double subregionHeight = (double)(imageHeight - 2 * yBoundary) / (double)ySubregionNum;

const double maxTrackDis = 100;
const int winSize = 21;

GoodFeaturesToTrackDetector_OCL oclFeatureDetector;
PyrLKOpticalFlow oclOpticalFlow;

vector<Point2f> *featuresCur = new vector<Point2f>();
vector<Point2f> *featuresLast = new vector<Point2f>();
vector<Point2f> featuresSub;
vector<unsigned char> featuresStatus;

int featuresIndFromStart = 0;
vector<int> featuresInd;

int totalFeatureNum = 0;
int subregionFeatureNum[2 * totalSubregionNum] = {0};

pcl::PointCloud<ImagePoint>::Ptr imagePointsCur(new pcl::PointCloud<ImagePoint>());
pcl::PointCloud<ImagePoint>::Ptr imagePointsLast(new pcl::PointCloud<ImagePoint>());

ros::Publisher *imagePointsLastPubPointer;
ros::Publisher *imageShowPubPointer;
cv_bridge::CvImage bridge;

static void download(const oclMat& ocl_mat, vector<Point2f>& vec)
{
  vec.resize(ocl_mat.cols);
  Mat mat(1, ocl_mat.cols, CV_32FC2, (void*)&vec[0]);
  ocl_mat.download(mat);
}

static void download(const oclMat& ocl_mat, vector<unsigned char>& vec)
{
  vec.resize(ocl_mat.cols);
  Mat mat(1, ocl_mat.cols, CV_8UC1, (void*)&vec[0]);
  ocl_mat.download(mat);
}

void imageDataHandler(const sensor_msgs::Image::ConstPtr& imageData) 
{
  timeLast = timeCur;
  timeCur = imageData->header.stamp.toSec() - 0.1163;

  cv_bridge::CvImageConstPtr imageDataCv = cv_bridge::toCvShare(imageData, "mono8");

  if (!systemInited) {
    remap(imageDataCv->image, image0, mapx, mapy, CV_INTER_LINEAR);
    oclImage0 = oclMat(image0);
    systemInited = true;

    return;
  }

  Mat imageLast, imageCur;
  oclMat oclImageLast, oclImageCur;
  if (isOddFrame) {
    remap(imageDataCv->image, image1, mapx, mapy, CV_INTER_LINEAR);
    oclImage1 = oclMat(image1);

    imageLast = image0;
    imageCur = image1;

    oclImageLast = oclImage0;
    oclImageCur = oclImage1;
  } else {
    remap(imageDataCv->image, image0, mapx, mapy, CV_INTER_LINEAR);
    oclImage0 = oclMat(image0);

    imageLast = image1;
    imageCur = image0;

    oclImageLast = oclImage1;
    oclImageCur = oclImage0;
  }
  isOddFrame = !isOddFrame;

  resize(oclImageLast, oclImageShow, showSize);
  oclImageShow.download(imageShow);
  cornerHarris(oclImageShow, oclHarrisLast, 2, 3, 0.04);
  oclHarrisLast.download(harrisLast);

  vector<Point2f> *featuresTemp = featuresLast;
  featuresLast = featuresCur;
  featuresCur = featuresTemp;

  pcl::PointCloud<ImagePoint>::Ptr imagePointsTemp = imagePointsLast;
  imagePointsLast = imagePointsCur;
  imagePointsCur = imagePointsTemp;
  imagePointsCur->clear();

  int recordFeatureNum = totalFeatureNum;
  detectionCount = (detectionCount + 1) % (detectionSkipNum + 1);
  if (detectionCount == detectionSkipNum) {
    oclMat oclFeaturesSub;
    for (int i = 0; i < ySubregionNum; i++) {
      for (int j = 0; j < xSubregionNum; j++) {
        int ind = xSubregionNum * i + j;
        int numToFind = maxFeatureNumPerSubregion - subregionFeatureNum[ind];

        if (numToFind > maxFeatureNumPerSubregion / 5.0) {
          int subregionLeft = xBoundary + (int)(subregionWidth * j);
          int subregionTop = yBoundary + (int)(subregionHeight * i);
          Rect subregion = Rect(subregionLeft, subregionTop, (int)subregionWidth, (int)subregionHeight);

          oclFeatureDetector.maxCorners = numToFind;
          oclFeatureDetector(oclImageLast(subregion), oclFeaturesSub);
          if (oclFeaturesSub.cols > 0) {
            oclFeatureDetector.downloadPoints(oclFeaturesSub, featuresSub);
            numToFind = featuresSub.size();
          } else {
            numToFind = 0;
          }

          int numFound = 0;
          for(int k = 0; k < numToFind; k++) {
            featuresSub[k].x += subregionLeft;
            featuresSub[k].y += subregionTop;

            int xInd = (featuresSub[k].x + 0.5) / showDSRate;
            int yInd = (featuresSub[k].y + 0.5) / showDSRate;

            if (harrisLast.at<float>(yInd, xInd) > 1e-7) {
              featuresLast->push_back(featuresSub[k]);
              featuresInd.push_back(featuresIndFromStart);

              numFound++;
              featuresIndFromStart++;
            }
          }
          totalFeatureNum += numFound;
          subregionFeatureNum[ind] += numFound;
        }
      }
    }
  }

  if (totalFeatureNum > 0) {
    Mat featuresLastMat(1, totalFeatureNum, CV_32FC2, (void*)&(*featuresLast)[0]);
    oclMat oclFeaturesLast(featuresLastMat);
    oclMat oclFeaturesCur, oclFeaturesStatus;

    oclOpticalFlow.sparse(oclImageLast, oclImageCur, oclFeaturesLast, oclFeaturesCur, oclFeaturesStatus);
    if (oclFeaturesCur.cols > 0 && oclFeaturesCur.cols == oclFeaturesStatus.cols) {
      download(oclFeaturesCur, *featuresCur);
      download(oclFeaturesStatus, featuresStatus);
      totalFeatureNum = featuresCur->size();
    } else {
      totalFeatureNum = 0;
    }
  }

  for (int i = 0; i < totalSubregionNum; i++) {
    subregionFeatureNum[i] = 0;
  }

  ImagePoint point;
  int featureCount = 0;
  for (int i = 0; i < totalFeatureNum; i++) {
    double trackDis = sqrt(((*featuresLast)[i].x - (*featuresCur)[i].x) 
                    * ((*featuresLast)[i].x - (*featuresCur)[i].x)
                    + ((*featuresLast)[i].y - (*featuresCur)[i].y) 
                    * ((*featuresLast)[i].y - (*featuresCur)[i].y));

    if (!(trackDis > maxTrackDis || (*featuresCur)[i].x < xBoundary || 
      (*featuresCur)[i].x > imageWidth - xBoundary || (*featuresCur)[i].y < yBoundary || 
      (*featuresCur)[i].y > imageHeight - yBoundary) && featuresStatus[i]) {

      int xInd = (int)(((*featuresLast)[i].x - xBoundary) / subregionWidth);
      int yInd = (int)(((*featuresLast)[i].y - yBoundary) / subregionHeight);
      int ind = xSubregionNum * yInd + xInd;

      if (subregionFeatureNum[ind] < maxFeatureNumPerSubregion) {
        (*featuresCur)[featureCount] = (*featuresCur)[i];
        (*featuresLast)[featureCount] = (*featuresLast)[i];
        featuresInd[featureCount] = featuresInd[i];

        point.u = -((*featuresCur)[featureCount].x - kImage[2]) / kImage[0];
        point.v = -((*featuresCur)[featureCount].y - kImage[5]) / kImage[4];
        point.ind = featuresInd[featureCount];
        imagePointsCur->push_back(point);

        if (i >= recordFeatureNum) {
          point.u = -((*featuresLast)[featureCount].x - kImage[2]) / kImage[0];
          point.v = -((*featuresLast)[featureCount].y - kImage[5]) / kImage[4];
          imagePointsLast->push_back(point);
        }

        featureCount++;
        subregionFeatureNum[ind]++;
      }
    }
  }
  totalFeatureNum = featureCount;
  featuresCur->resize(totalFeatureNum);
  featuresLast->resize(totalFeatureNum);
  featuresInd.resize(totalFeatureNum);

  sensor_msgs::PointCloud2 imagePointsLast2;
  pcl::toROSMsg(*imagePointsLast, imagePointsLast2);
  imagePointsLast2.header.stamp = ros::Time().fromSec(timeLast);
  imagePointsLastPubPointer->publish(imagePointsLast2);

  showCount = (showCount + 1) % (showSkipNum + 1);
  if (showCount == showSkipNum) {
    bridge.image = imageShow;
    bridge.encoding = "mono8";
    sensor_msgs::Image::Ptr imageShowPointer = bridge.toImageMsg();
    imageShowPubPointer->publish(imageShowPointer);
  }
}

int main(int argc, char* argv[])
{
  ros::init(argc, argv, "featureTracking_ocl");
  ros::NodeHandle nh;

  Size imageSize = Size(imageWidth, imageHeight);
  mapx.create(imageSize, CV_32FC1);
  mapy.create(imageSize, CV_32FC1);
  initUndistortRectifyMap(kMat, dMat, Mat(), kMat, imageSize, CV_32FC1, mapx, mapy);

  oclOpticalFlow.winSize = Size(winSize, winSize);

  ros::Subscriber imageDataSub = nh.subscribe<sensor_msgs::Image>("/image/raw", 1, imageDataHandler);

  ros::Publisher imagePointsLastPub = nh.advertise<sensor_msgs::PointCloud2> ("/image_points_last", 5);
  imagePointsLastPubPointer = &imagePointsLastPub;

  ros::Publisher imageShowPub = nh.advertise<sensor_msgs::Image>("/image/show", 1);
  imageShowPubPointer = &imageShowPub;

  ros::spin();

  return 0;
}
