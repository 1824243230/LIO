
#include <csignal>
#include <exception>
#include <glog/logging.h>
#include <ros/ros.h>
#include "lio/super_lio_reloc.h"
#include "ros/ROSWrapper.h"


using namespace LI2Sup;

void SigHandle(int sig) {
  g_flag_run = false;
}

int main(int argc, char** argv){
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = true;
  ros::init(argc, argv, "lio");
  signal(SIGINT, SigHandle);
  ros::NodeHandle nh;
  try {
    LoadParamFromRos(nh);
  } catch (const std::exception& error) {
    ROS_FATAL_STREAM("LIO configuration error: " << error.what());
    google::ShutdownGoogleLogging();
    return 1;
  }

  ROSWrapper::Ptr data_wrapper = std::make_shared<ROSWrapper>();
  auto lio = std::make_shared<SuperLIOReLoc>();
  lio->setROSWrapper(data_wrapper);
  try {
    lio->init();
  } catch (const std::exception& error) {
    ROS_FATAL_STREAM("LIO initialization failed: " << error.what());
    google::ShutdownGoogleLogging();
    return 1;
  }

  ros::Rate rate(500);  // 500 Hz
  while (ros::ok() && g_flag_run) {
    data_wrapper->spinOnce();
    lio->process();
    rate.sleep();
  }

  lio->saveMap();
  lio->printTimeRecord();
  google::ShutdownGoogleLogging();
  return 0;
}