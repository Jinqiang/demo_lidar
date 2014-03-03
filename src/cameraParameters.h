#ifndef camera_parameters_h
#define camera_parameters_h

#include <opencv/cv.h>
#include <sensor_msgs/Image.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/features2d/features2d.hpp>
#include "opencv2/contrib/contrib.hpp"
#include "opencv2/calib3d/calib3d.hpp"

#include <cv_bridge/cv_bridge.h>

const int imageWidth = 744;
const int imageHeight = 480;

double kImage[9] = {4.177343016733e+002, 0, 3.715643918956e+002, 
                    0, 4.177970397634e+002, 1.960688121183e+002, 
                    0, 0, 1};

double dImage[4] = {-3.396867101163e-001, 1.309347902588e-001, -2.346791258754e-004, 2.209387016957e-004};

#endif

