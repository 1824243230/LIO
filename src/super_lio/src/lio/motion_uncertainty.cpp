#include "lio/motion_uncertainty.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace LI2Sup::motion {
M3 hat(const V3& v) {
  M3 K; K << 0,-v.z(),v.y(),v.z(),0,-v.x(),-v.y(),v.x(),0; return K;
}
void integrals(const V3& w, double dt, M3& velocity, M3& position) {
  const V3 phi = w*dt;
  const double t2 = phi.squaredNorm();
  double a, b, c;
  // 小角度时闭式分子接近零，用级数避免相减损失精度。
  if (t2 < 1e-4) {
    a = 0.5-t2/24+t2*t2/720;
    b = 1.0/6-t2/120+t2*t2/5040;
    c = 1.0/24-t2/720+t2*t2/40320;
  } else {
    const double t = std::sqrt(t2);
    a = (1-std::cos(t))/t2;
    b = (t-std::sin(t))/(t2*t);
    c = (0.5*t2+std::cos(t)-1)/(t2*t2);
  }
  const M3 K = hat(phi), K2 = K*K;
  velocity = dt*(M3::Identity()+a*K+b*K2);
  position = dt*dt*(0.5*M3::Identity()+b*K+c*K2);
}
void relativeCovariances(std::vector<Knot>& knots) {
  if (knots.empty()) return;
  const auto& end = knots.back();
  // 反向循环开始时：δx_end = phi * δx_i + η，noise = Cov(η)。
  // 每次先计算当前节点相对末端的不确定性，再把递推范围扩展到前一节点。
  M18 phi = M18::Identity(), noise = M18::Zero();
  for (size_t i = knots.size(); i-- > 0;) {
    auto& k = knots[i];
    const V3 t = end.R.transpose()*(k.p-end.p);
    Eigen::Matrix<double,6,18> A = Eigen::Matrix<double,6,18>::Zero(), B = A;
    A.block<3,3>(0,0) = end.R.transpose()*k.R;
    A.block<3,3>(3,3) = end.R.transpose();
    B.block<3,3>(0,0) = -M3::Identity();
    B.block<3,3>(3,0) = hat(t);
    B.block<3,3>(3,3) = -end.R.transpose();
    // 同一先验引起的误差必须先相加再求协方差；独立相加 Pi/Pe 会漏掉互相关。
    const Eigen::Matrix<double,6,18> joint = A+B*phi;
    k.relative = joint*k.covariance*joint.transpose()+B*noise*B.transpose();
    k.relative = (0.5*(k.relative+k.relative.transpose())).eval();
    // Same random variable at scan end: make exact cancellation explicit.
    if (i == knots.size()-1) k.relative.setZero();
    noise += phi*k.noise*phi.transpose();
    phi = (phi*k.transition).eval();
  }
}
M3 pointCovariance(const V3& point, double time, const std::vector<Knot>& knots) {
  // NaN 与所有时间比较均为 false，会使 upper_bound 返回 end()。
  // 非法查询返回非有限协方差，由 Observe 拒绝，不能伪装为零噪声。
  if (!std::isfinite(time) || !point.allFinite())
    return M3::Constant(std::numeric_limits<double>::quiet_NaN());
  if (knots.size()<2 || time>=knots.back().time) return M3::Zero();
  auto next = std::upper_bound(knots.begin(),knots.end(),time,
      [](double t,const Knot& k){return t<k.time;});
  const Knot* head = &knots.front();
  double u = 0;
  if (next != knots.begin()) {
    head = &*std::prev(next);
    u = std::clamp((time-head->time)/(next->time-head->time),0.0,1.0);
  } else next = knots.begin();
  // 两端协方差已表达在同一扫描末端坐标系，可做保持 PSD 的凸插值。
  const M6 covariance = (1-u)*head->relative+u*next->relative;
  const V3 translation = knots.back().R.transpose()*
      ((1-u)*head->p+u*next->p-knots.back().p);
  Eigen::Matrix<double,3,6> J;
  J.leftCols<3>() = -hat(point-translation);
  J.rightCols<3>() = M3::Identity();
  const M3 result = J*covariance*J.transpose();
  return 0.5*(result+result.transpose());
}
}
