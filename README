Depth Enhanced Monocular Odometry (Demo) is a monocular visual odometry method assisted by depth maps, and optionally an IMU. The program contains three major nodes. The "featureTracking" node extracts and tracks Harris corners by Kanade Lucas Tomasi (KLT) tracker provided in the OpenCV library. The "visualOdometry" node takes the tracked features and estimates frame to frame motion. Features associated with depth (either from the depth map or triangulated from previously estimated camera motion) are used to solve the 6DOF motion, and features without depth help solve orientation. The "bundleAdjust" node refines the estimated motion with bundle adjustment. It processes sequences of images using Incremental Smoothing and Mapping (iSAM) open source library. 

The program is tested on ROS Indigo, on a laptop computer with 2.5 GHz quad cores and 6 Gib memory. This version uses a camera with point cloud perceived by a 3D lidar.

Wiki Webpage: http://wiki.ros.org/demo_lidar

Another version of the program that uses an RGBD camera is available at

Wiki Webpage: http://wiki.ros.org/demo_rgbd

GitHub Code: https://github.com/jizhang-cmu/demo_rgbd.git

Also, a slimmed down version (without the bundle adjustment) that uses an RGBD camera is available for running on embedded systems

GitHub Code: https://github.com/jizhang-cmu/demo_rgbd_slimmed.git

How to use:

1) Install additional packages by "sudo apt-get install libsuitesparse-dev libeigen3-dev libsdl1.2-dev". These packages are required by iSAM library for the bundle adjustment.

2) Download the program file to a ROS directory, unpack the file and rename the folder to “demo_lidar” (GitHub may add "-xxx" to the end of the folder name). Use “catkin_make” to build the program, then “roslaunch demo_lidar.launch”. The launch file should start the program and rviz.

3) Download datasets from the following website. Make sure the data files are for the camera and lidar (not RGBD camera) version. Play the data files with “rosbag play data_file_name.bag”. Note that if a slow computer is used, users can try to play the data files at a lower speed, e.g. “rosbag play data_file_name.bag -r 0.5” plays the data file at half speed.

Datasets can be downloaded at: http://www.frc.ri.cmu.edu/~jizhang03/projects.htm

Notes: 
  
Camera intrinsic parameters (K and D matrices) are defined in the "src/cameraParameters.h" file.

It is possible to accelerate feature tracking with a GPU. To do this, simply replace the "src/featureTracking.cpp" file with the "src/featureTracking_ocl.cpp" file and recompile. We use OpenCL to communicate with the GPU. Users first need to install OpenCL in a version no earlier than 1.1 with full profile, then install OpenCV using CMake with flag WITH_OPENCL=ON.
