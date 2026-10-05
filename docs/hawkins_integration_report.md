# Hawkins Integration Report — 阶段 4

验证日期：2026-10-05。对象为 SubT-MRS Hawkins Long Corridor RC 的原版 Super-LIO；退化检测等增强项全部关闭。完整 279.997 s bag 已原速回放两次，第二次使用不附加外参参数的 `Hawkins.launch` 默认入口。

## 1. 阅读的文件与来源

- 本工程：`src/super_lio/launch/Hawkins.launch`、`src/super_lio/config/Hawkins.yaml`、`src/super_lio/src/ros/ROSWrapper.cpp`、`src/super_lio/src/lio/super_lio.cpp`、`src/basic/include/basic/alias.h`、`src/basic/src/Manifold.cpp`、`src/super_lio/rviz/lio.rviz`。
- 本地数据：`/home/hyh/lidarbag/Long_Corridor_Rosbag/raw_data_core_2023-07-25-03-01-44.bag` 和 `ground_truth_path.csv`。
- [挑战赛硬件说明](https://github.com/yuanjun-gao/ICCV2023_SLAM_Challenge/blob/master/Hardware_Information.md)：雷达为 Velodyne VLP-16（10 Hz），IMU 为 Epson M-G365（200 Hz）。[Long Corridor RC 官方页面](https://superodometry.com/iccv23_challenge_LiI)列出专属[外参文件](https://drive.google.com/file/d/1bB3jfEJeTf_XoLUHKOaxCNF_MCkiQTol/view)及[内参文件](https://drive.google.com/file/d/10rv5dg5un7kUveTPS3XBx8IuIYGg9r2h/view)。取得的“内参”文件仅包含 RGB 相机模型，不包含 VLP-16 逐线修正。

## 2. 数据流与关键接口

`/velodyne_packets` (`velodyne_msgs/VelodyneScan`) → Hawkins launch 内仅启用的 `velodyne_pointcloud/TransformNodelet` → `/velodyne_points` (`sensor_msgs/PointCloud2`) → `ROSWrapper::stdMsgHandler` 的 `VELO16` 分支 → `finalizeLidarScan` → 去畸变、降采样、观测、ESKF 更新 → OctVox 地图。`/imu/data` (`sensor_msgs/Imu`) 直接进入 `ROSWrapper::imuHandler`。

Hawkins launch 的解码器使用 `model=VLP16` 与 ROS 安装包的通用 `VLP16db.yaml`。不会启动实时 UDP driver；如果单独启动其他解码器，可设 `decode_packets:=false` 避免重复发布。LiDAR/IMU topic 分别为 `/velodyne_points`、`/imu/data`，`lidar_type=3`，雷达过滤率 3、盲区 0.4 m、最大距离 100 m、体素大小 0.5 m。IMU 噪声仍是标明为初始调参值的配置，并非官方 Allan 标定。

## 3. 外参及坐标约定

官方 YAML 的键为 `laser_to_imu`，平移为 `[0.08, 0.029, 0.03]` m。代码在去畸变前计算 `R * p_L + t`，因此实际需要的是 `p_I = R_IL p_L + t_IL`；参数数组先放平移，再放旋转的 Eigen **列主序**九项。没有使用 M2DGR 或单位外参。

官方文件的原始旋转矩阵是：

```text
[ 0.999212900, -0.000519121,  0.004000000 ]
[ 0.000516111,  0.999218492, -0.000939132 ]
[-0.004000000,  0.000802565,  0.999993652 ]
```

其行列式为 `0.998442677`，`||RᵀR-I||_F=0.002213966`，不能直接作为严格的 SO(3) 旋转；本工程 `SE3::inverse()` 使用转置作为逆。为保证流形运算，`Hawkins.launch` 默认使用该官方矩阵的最近正交旋转：`R = U diag(1,1,det(UVᵀ)) Vᵀ`，其中 `UΣVᵀ` 是原始矩阵的 SVD。投影改变量 `||R_projected-R_official||_F=0.001107417`。完整转换后的 12 个数值写在专用 launch 的 `lidar_imu` 参数中，可在命令行覆盖。**此投影是对官方文件数值的数学修正，仍应向数据集方确认原始矩阵为何不严格正交。**

## 4. 时间戳与字段

bag 录制时间为 2023 年，传感器消息 `header.stamp` 为 2018 年；Super-LIO 根据消息时间同步，因此 `Hawkins.launch` 设 `use_sim_time=false`。回放的 `--clock` 不参与 LIO 时间计算。

实测解码字段为 `x/y/z/intensity/time: FLOAT32`、`ring: UINT16`；无 `t`、`timestamp`。`time` 是从扫描起点计的相对秒数，抽样范围 `0–0.10083945 s`。点云起点取 `PointCloud2.header.stamp`，扫描终点取起点加有效点最大相对时间，去畸变查询时间为 `scan.start_time + point.offset_time`。两次回放未发现点云或 IMU header 时间倒退；逐点时间抽样无非有限值。Hawkins 输入层会拒绝缺失 `FLOAT32 time` 字段的点云。

## 5. 实际运行结果

在独立 ROS master（端口 11412）上，第一轮以官方外参投影值作为临时命令行参数运行；第二轮直接执行默认 `roslaunch super_lio Hawkins.launch`。使用以下命令回放完整时间段，过滤掉不参与 LIO 的相机话题：

```bash
source /opt/ros/noetic/setup.bash
source /home/hyh/rong_ws/Super-LIO/devel/setup.bash
roslaunch super_lio Hawkins.launch
# 另一终端：
rosbag play /home/hyh/lidarbag/Long_Corridor_Rosbag/raw_data_core_2023-07-25-03-01-44.bag --clock -q --topics /velodyne_packets /imu/data
```

两轮均完成全长回放；未更改核心算法。下面是第二轮的独立订阅统计：

| 话题 | 收到消息 | 按消息时间计算的平均频率 | 时间倒退 |
| --- | ---: | ---: | ---: |
| `/velodyne_packets` | 2777 | 9.915 Hz | 0 |
| `/imu/data` | 55957 | 199.841 Hz | 0 |
| `/velodyne_points` | 2776 | 9.915 Hz | 0 |
| `/lio/odom` | 2769 | 9.915 Hz | 0 |
| `/lio/cloud_world` | 553 | 1.983 Hz | 0 |

第二轮的统计订阅器错过最初一帧解码点云；第一轮完整收到 2777/2777 帧。启动日志确认 `KF init done`、`Map init done`，之后 `/lio/odom` 持续至 bag 末端，位姿非有限值计数为 0，最大相邻里程计平移约 0.376 m。RViz 显示连贯的走廊点云与轨迹，见 [截图](hawkins_stage4_rviz.png)。`/lio/cloud_world` 持续输出；结合 `stateProcess` 代码路径与没有观测拒绝日志，可以推断正常观测后的 OctVox 地图插入持续进行。本次未增加 OctVox 内部容量遥测，无法直接报告体素数量。

第二轮从 `/lio/path` 导出 2229 个轨迹点（该话题仅在移动超过 0.1 m 时加入路径）。与本地真值按 0.02 s 对齐，共匹配 1111 对；`evo_ape tum -r trans_part --align` 使用 SE(3) 对齐、**不拟合尺度**，ATE RMSE 为 **1.503 m**，均值 1.198 m，最大值 5.922 m。这支持“未见明显发散”的判断，但不是对外参和噪声模型精度的最终认证。轨迹、转换后的真值、统计和 APE 结果保存在 `/home/hyh/rong_ws/results/hawkins_stage4/`，复算命令为：

```bash
evo_ape tum /home/hyh/rong_ws/results/hawkins_stage4/groundtruth.tum /home/hyh/rong_ws/results/hawkins_stage4/estimate.tum -r trans_part --align --t_max_diff 0.02
```

## 6. 错误、修改、影响与编译

- 解码器偶发 `Packet containing angle overflow`，角度从约 359.7° 跨过 0°；两轮点云、里程计持续输出，无对应缺帧或时间倒退。启动时还有 `No Azimuth Cache configured for model VLP16` 提示。未见 Super-LIO 的时间同步错误、点时间异常、观测拒绝日志或 NaN。
- 阶段 4 修改 `Hawkins.launch` 的专属外参默认值及来源说明、`Hawkins.yaml` 的外参说明；新增本报告和 RViz 截图。**M2DGR.launch、其他数据集配置、ESKF/IESKF、IMU propagation、OctVox、HKNN、measurement 均未修改。**
- 本阶段的增量 `catkin_make -j2 -l2` 通过；阶段 3 的 13/13 项测试通过。本阶段仅改配置和文档，`roslaunch --files`、`--dump-params`、`git diff --check` 通过。默认 launch 已经完成实际全长运行。

## 7. 当前遗留问题与下一阶段建议

1. 官方提供的 LiDAR 外参旋转并非严格正交；目前使用最近 SO(3) 投影，建议向数据集方核对原始值和坐标系定义。
2. 官方内参链接只提供 RGB 相机参数，当前 VLP-16 使用 ROS 通用逐线校正；若取得本机专属 Velodyne calibration，应在 Hawkins launch 中替换 `velodyne_calibration`。
3. IMU 白噪声与偏置随机游走仍是初始值，需用该 Epson M-G365 的正式 Allan 参数做后续精度评估。
4. 若需长期地图容量、观测成功率及资源占用统计，可在独立评测阶段增加非侵入式监测；本阶段没有修改 OctVox 或估计器。

阶段 4 到此停止，不自动进入退化增强阶段。
