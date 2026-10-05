# Hawkins baseline 与几何增强入口

`Hawkins.launch` 保持原版 Super-LIO：Hawkins 数据解码、话题、外参及 IMU 参数不变，五项几何增强默认关闭。

`degeneracy_Hawkins.launch` 包含 `Hawkins.launch`，再加载现有 `config/geometry.yaml` 并覆盖开关。默认开启退化检测、谱可靠性、bump layer、informed sampling 和 bump measurement；运动传播与不确定性运动补偿仍按基础入口关闭。该入口没有复制 M2DGR 外参或更改状态估计、地图、HKNN、IMU 传播代码。

```text
Hawkins.launch ───────────────────────────────→ 原版 Super-LIO (baseline)
        │
        └─ degeneracy_Hawkins.launch
             + geometry.yaml
             + 五项几何开关 ────────────────────→ 几何增强 Super-LIO
```

## 运行

```bash
source /opt/ros/noetic/setup.bash
source /home/hyh/rong_ws/Super-LIO/devel/setup.bash
roslaunch super_lio Hawkins.launch
# 或（另一轮实验）
roslaunch super_lio degeneracy_Hawkins.launch
```

另一个终端回放：

```bash
rosbag play /home/hyh/lidarbag/Long_Corridor_Rosbag/raw_data_core_2023-07-25-03-01-44.bag --clock --topics /velodyne_packets /imu/data
```

两个入口不能在同一个 ROS master 上同时启动，因为它们使用相同节点名与输出话题。做对照实验时分别启动，并保持同一 bag、外参与参数。`Hawkins.launch` 内已经包含 packet→PointCloud2 解码器；如果另行启动解码器，两种入口都可用 `decode_packets:=false`。`rviz:=false` 可用于无界面评测，`velodyne_calibration:=...` 可替换逐线雷达校正文件。

增强入口的 `enable:=true` 控制退化检测、bump layer、informed sampling、bump measurement；`spectral:=true` 控制谱可靠性。设置 `enable:=false spectral:=true` 时只启用退化检测与谱可靠性；两个都为 `false` 时关闭五项增强。谱可靠性需要退化检测，因此 launch 在 `spectral:=true` 时仍会开启检测。`csv_path:=/absolute/path/geometry.csv` 可记录逐帧诊断；目标目录须预先存在，默认不写 CSV。

诊断话题位于 `/super_lio/`，包括 `degeneracy_status`、`weak_dim`、`spectral_valid`、`spectral_attenuated_modes`、`bump_candidate_count`、`bump_accepted_count`、`bump_fallback` 等。bump sampling/measurement 只会在历史地图中找到支持的局部表面且通过质量门控时产生约束；开启开关不保证每帧都有有效 bump 观测。

使用的几何阈值来自现有 `geometry.yaml`，尚未针对 Hawkins 调参；因此“增强入口”表示五项功能被启用，不预先声称相对 baseline 的精度必然提升。Hawkins 原始传感器和外参限制见 [集成报告](hawkins_integration_report.md)。

## 验证

`roslaunch --files` 和 `--dump-params` 显示 baseline 的五项开关均为 `false`，增强入口均为 `true`，其他 Hawkins 传感器参数一致。30 秒实际回放完成 IMU/KF 和地图初始化，产生 290 行诊断：退化状态变化，`spectral_valid=1`，最多衰减 3 个模态，bump 候选最多 1031、接受最多 200，说明增强路径在运行。

完整 279 s 序列也完成回放，日志确认 `KF init done` 和 `Map init done`，未见 LIO 同步、NaN 或观测拒绝错误。记录了 2769 行诊断，其中退化状态为正常/弱/退化的帧数分别是 11/1467/1291；`spectral_valid=1` 共 2769 帧，至少衰减一个模态的有 2765 帧；有 bump 候选的 2758 帧，接受 bump 约束的 2745 帧，13 帧使用回退。诊断数据无非有限值。详见 `/home/hyh/rong_ws/results/hawkins_degeneracy/geometry.csv`。

使用与 [baseline 集成报告](hawkins_integration_report.md)相同的真值、`evo_ape tum -r trans_part --align --t_max_diff 0.02`，且不拟合尺度：

| 入口 | 匹配的真值对数 | ATE RMSE | ATE 最大值 |
| --- | ---: | ---: | ---: |
| `Hawkins.launch` baseline | 1111 | 1.503 m | 5.922 m |
| `degeneracy_Hawkins.launch` 五项全开 | 1105 | 1.860 m | 8.102 m |

**本次配置在 Hawkins Long Corridor 上未取得精度提升**；RMSE 约高 0.357 m（23.8%）。两个入口均可运行，增强指标也证实五项功能进入了实际处理路径，但 `geometry.yaml` 是通用阈值，尚未针对该场景标定。后续如需改善精度，应在固定相同输入与外参后做逐项消融和阈值评估；不能把“功能开启”当作“精度改善”。增强版轨迹与 APE 结果保存在 `/home/hyh/rong_ws/results/hawkins_degeneracy/`。本次只新增 launch 与文档，`catkin_make -j2 -l2` 通过，未修改原版 Hawkins/M2DGR 配置或任何估计核心。
