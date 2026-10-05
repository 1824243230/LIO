# Hawkins Long Corridor 适配：阶段 1 核查

范围：只核对原始 bag、已安装 ROS1 解码器和 Super-LIO 输入接口。本阶段没有改动工程配置或状态估计算法，也没有启动数据集回放。

## 已核实的数据与接口

- 本地 bag：`/home/hyh/lidarbag/Long_Corridor_Rosbag/raw_data_core_2023-07-25-03-01-44.bag`；`rosbag info` 显示约 279 s，`/velodyne_packets` 2777 帧、`/imu/data` 55957 帧、`/camera_1/image_raw` 6720 帧。本次 LIO 适配只需要前两种消息。
- 主办方的 [硬件说明](https://github.com/yuanjun-gao/ICCV2023_SLAM_Challenge/blob/master/Hardware_Information.md)将该载荷的雷达列为 Velodyne VLP-16（10 Hz），IMU 列为 Epson M-G365（200 Hz）。抽样 bag 消息的 frame 分别为 `cmu_rc2_velodyne` 和 `epson`；首帧雷达包含 76 个 1206 字节 packet，packet 产品标识为 `0x22`。
- 主办方的 [ROS1 转换说明](https://github.com/yuanjun-gao/ICCV2023_SLAM_Challenge/blob/master/README.md#instructions-for-running-velodyne-driver)推荐使用 `velodyne_pointcloud` 将 `VelodyneScan` 转为 `PointCloud2`。本机 ROS Noetic 已安装 `velodyne_msgs`、`velodyne_pointcloud` 及 VLP16 解码参数。安装的 `PointcloudXYZIRT` 输出 `x/y/z/intensity/ring/time`；Super-LIO 的 `VELO16` 输入已有 `x/y/z/intensity/time` 路径，时间单位直接按秒使用，最终经 `finalizeLidarScan` 取真实点时间边界。
- 原始 bag 只有 `/velodyne_packets`，没有 `/velodyne_points`；仅改雷达 topic 无法使原工程直接订阅。后续独立 `Hawkins.launch` 应在 bag 回放时只加载 packet→pointcloud 解码 nodelet，再让原 Super-LIO 订阅 `/velodyne_points` 和 `/imu/data`。现成的 `VLP16_points.launch` 还会启动实时 UDP driver 和 LaserScan nodelet，不宜原样作为仅回放 bag 的入口。

## 时间与标定待核实

- bag 录制/播放时间约为 `1690254112`（2023），而首帧 IMU/LiDAR `header.stamp` 约为 `1517157219`（2018）。在 bag 开始、约 100 s 和约 270 s 的抽样中，这个约 `173096894 s` 的偏移保持稳定。IMU 与雷达**消息头时间**位于同一时间轴；扫描首包与末包相差约 0.1 s。`rosbag play --clock` 发布的是 bag 时间，不能把它当成传感器头时间。阶段 2 须实测解码后点云时间字段、Super-LIO 同步及 RViz 展示；不要用错误的全局时间偏移去改滤波器。
- 本地数据目录只有 bag、`ground_truth_path.csv` 和地图图片，没有外参 YAML。[主办方专门提供了 Long Corridor RC 外参文件](https://drive.google.com/file/d/1bB3jfEJeTf_XoLUHKOaxCNF_MCkiQTol/view?usp=drive_link)，但当前网页只显示 Google Drive 登录页，尚未取得数值。**TODO：取得该官方文件，确认其中 LiDAR→IMU 变换的方向、坐标轴、单位及数组顺序，再填入 `lio/extrinsic/lidar_imu`。不能沿用 M2DGR 数值或假设单位外参。**
- 原工程 `LoadParamFromRos` 要求 `lidar_imu` 恰为 12 个有限数（平移 3 项、旋转 9 项），内部构造 `SE3(R,t)`；运行前必须明确官方标定的矩阵储存顺序。重力模长及 IMU 噪声也应采用 Hawkins 官方参数或经过验证的实测值，不能直接认定 M2DGR 参数正确。

## 下一阶段边界

取得/核实官方 LiDAR–IMU 标定后，新增独立 `Hawkins.yaml` 与 `Hawkins.launch`，以 bag 短片验证解码点字段、时间顺序、IMU 同步和原版里程计输出。`M2DGR.launch` 及 ESKF、IMU propagation、OctVox、HKNN 和 point-to-plane measurement 保持原样。退化增强入口 `degeneracy_Hawkins.launch` 留待基础入口验证后再实施。
