# SE(3)-LIO 思路在 Super-LIO 中的适配

## 审查结论

参考本地 `/home/hyh/LIO/se3-lio/cpp/se3_lio/core/state_predict.cpp` 的
`predictStateWithCov`、`calculateRelPoseWithCov`、`calculateUndistCloudCov`，以及
[SE(3)-LIO 原论文](https://arxiv.org/abs/2603.16118)。

修改前的 Super-LIO 使用 SO(3) 姿态更新、世界系常加速度位置更新、扫描内姿态 slerp 和二次位置插值。
虽然存在 SO(3)/SE(3) 类型，但不等于已采用 SE(3)-LIO 的创新。传播只保留当前状态协方差，
去畸变没有位姿联合分布，平面观测使用固定 0.001 m² 方差。

| 机制 | 本次处理 |
|---|---|
| 原论文 SE(3) 指数传播、机体系速度、17 维误差 | 不直接移植。与本项目世界系速度、右旋转扰动、世界平移扰动及 18 维重力状态不兼容 |
| 传播中考虑旋转对平移的影响 | 实现恒定机体系角速度/比力下的旋转积分，传播和扫描内插值使用同一公式 |
| jointly distributed poses | 在现有误差坐标下实现等价的一阶联合协方差递推，包含共享先验及跨时刻相关性 |
| uncertainty-aware motion compensation | 相对位姿协方差传播到去畸变点，再投影到平面/曲面残差方差 |

这是借鉴并适配，不是完整 SE(3)-LIO 复现，也不是新的首创声明。未复制原仓库实现代码。

## 传播模型与坐标约定

保留 `δx=[δθ_body, δp_world, δv_world, δbg, δba, δg_world]`。
每个 IMU 区间沿用两端测量平均及扫描结束插值；令 `w` 为去偏角速度、`a` 为去偏比力：

```
R(t+h) = R(t) Exp(w h)
Iv(w,h) = integral_0^h Exp(w s) ds
Ip(w,h) = integral_0^h (h-s) Exp(w s) ds
p(t+h) = p(t) + v(t) h + R(t) Ip a + 0.5 g h²
v(t+h) = v(t) + R(t) Iv a + g h
```

Iv/Ip 使用闭式 SO(3) 积分及小角度级数。这里连续旋转的比力参与平移积分，
没有将世界系速度直接当作 SE(3) 的机体系 twist。姿态与位置连续，区间内速度为位置导数；
不保证 IMU 区间边界处加速度连续，也不称作全局平滑器。

`F` 和 `G` 同步增加位置对旋转、加速度偏置、重力及噪声的导数，
速度/位置对陀螺偏置的导数通过双精度中心差分积分矩阵获得（步长 1e-5 rad/s）。
保留原项目 `G Q Gᵀ` 的离散噪声参数约定，不改变噪声单位，不声称完成 IMU 噪声标定。
高频 IMU 输出也应用相同名义传播模型。

## 联合相对位姿协方差

扫描内只保存成功传播的严格递增节点。对节点 i 到扫描结束 e：

```
δxe = Φei δxi + ηei,  Cov(ηei)=Sei
Cov(δxi,δxe) = Pi Φeiᵀ
```

相对变换为 `Ri_rel=Reᵀ Ri, ti=Reᵀ(pi-pe)`。采用扫描末端坐标系中的
相对旋转左扰动及加性平移扰动：

```
δθrel = Reᵀ Ri δθi - δθe
δtrel = Reᵀ δpi + hat(ti) δθe - Reᵀ δpe
δξrel = Ai δxi + Bi δxe
Σrel = (Ai+Bi Φei) Pi (Ai+Bi Φei)ᵀ + Bi Sei Biᵀ
```

该写法显式保留相关性，并避免“大协方差相减”引起抵消误差。
从末端向前累计 Φ 和 S，节点计算复杂度 O(N)，不构造随扫描长度增长的稠密联合矩阵。
共享绝对平移误差抵消，但速度、偏置等导致的相对运动不确定性仍保留。
扫描末端相对自身的协方差显式设为零。

节点间 Σrel 采用同一末端坐标系中的凸插值，平移也线性插值用于点 Jacobian。
这是 IMU 间隔内的近似，不是每个点时刻重新传播得到的精确联合分布。
它保持半正定和末端零方差，但不保证相对精确连续模型保守；低 IMU 频率、极快运动需另行验证。

## 点云及观测融合

对于已经去畸变的点 q：

```
Jq = [-hat(q-ti), I]
Cq = Jq Σrel Jqᵀ
Rp = 0.001 + (Rcurrentᵀ n)ᵀ Cq (Rcurrentᵀ n)
```

原始 LiDAR 到 IMU 外参先作用到点，杆臂包含在 q 中；没有把当前绝对位姿协方差重复加到点噪声。
Cq 固定在本帧传播先验的末端 IMU 坐标系，观测迭代时用当前姿态旋转法向。
曲面残差采用其平移 Jacobian（不要求单位法向）做相同投影，追加到原曲面方差。
平面/曲面互斥替换减去每个点实际的平面权重；谱可靠性仍作用于最终观测系统。

当前中心体素采样和 informed sampling 都复制实际点。使用完整点云 XYZ 精确键关联 Cq，
不改 intensity、不借用 PCL padding。坐标完全重合的多个点使用协方差之和作为保守上界。
若以后换成质心体素滤波，必须同时实现协方差传递；缺失或非有限关联会拒绝本帧更新。
原采样回退重新取对应 Cq，不会悄悄回到固定权重。

局限：这里只保留逐点边缘协方差用于加权；点与点、点与滤波先验、点与历史地图之间仍有相关性，
没有实现完整相关测量滤波，不能宣称统计一致性。地图维护仍使用原有点地图，没有新增地图点协方差。
完全不可观方向仍需外部信息。新增全点协方差和关联表有内存及运算成本，未测量真实耗时。

## 开关与启动

所有普通配置默认关闭，现有 geometry 及 spectral 开关不变：

- `/lio/motion/smooth_propagation`：连续旋转比力传播。
- `/lio/motion/uncertainty_aware`：联合位姿去畸变不确定性。

M2DGR 显式入口（四种 smooth/uamc 组合均可使用）：

```bash
source /opt/ros/noetic/setup.bash
source /home/hyh/rong_ws/Super-LIO/devel/setup.bash
roslaunch super_lio degeneracy_M2DGR.launch enable:=false spectral:=false smooth:=true uamc:=true
```

仅 UAMC 时沿用基础名义传播及其 F/G（后续代码审查已补齐位置导数和位置噪声项）；仅 smooth 时不计算点协方差。
其他传感器启动文件可在启动前设置以上两个 ROS 参数。参数在节点初始化读取，修改后需重启。
为避免 ROS master 上的残留参数影响对照，使用显式入口传入 `smooth:=false uamc:=false`。

## 验证

新增 `motion_uncertainty` 测试：独立数值积分对比、零角速率极限、传播 Jacobian 有限差分、
扫描末端连续性与分段积分、共享误差抵消、速度/未来噪声贡献、完整联合矩阵对照、
协方差半正定、远距离旋转误差放大、非统一平面权重替换、实际 Observe 后验及原采样回退。

未运行数据集，不报告本次修改的 ATE/RPE、精度提升或实时性结论。

本轮结果：完整 `catkin_make -j2 -l2` 编译成功；9 项 CTest 全部通过，包含新增原始点云→传播→降采样关联检查；四种 smooth/uamc 组合通过 `roslaunch --dump-params` 检查（未启动节点）；`git diff --check` 通过。构建仍提示原环境 PCL 1.10/1.12 混链，本轮未改动依赖。

复查命令：

```bash
source /opt/ros/noetic/setup.bash
source /home/hyh/rong_ws/Super-LIO/devel/setup.bash
catkin_make -C /home/hyh/rong_ws/Super-LIO -j2 -l2
ctest --test-dir /home/hyh/rong_ws/Super-LIO/build --output-on-failure
```
