#include "lio/geometry.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace LI2Sup::geometry;
void check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
int main() {
  Options o; Diagnostics stats;
  M6 A=M6::Identity(); V6 b=V6::Ones();
  check(applySpectralReliability(A,b,o,stats) && A==M6::Identity() && b==V6::Ones(),
        "well-conditioned system is exactly unchanged");
  A(5,5)=0; b[5]=10;
  check(applySpectralReliability(A,b,o,stats) && A(5,5)==0 && b[5]==0,
        "nullspace has neither information nor residual forcing");
  check(stats.spectral_attenuated==1 && stats.spectral_min_gain==0, "null diagnostics");
  // Independent orthogonal chart with mixed rotation/translation eigenvectors.
  M6 U=M6::Identity();
  for(int k=0;k<5;++k) {
    M6 R=M6::Identity(); const double c=std::cos(0.3+k), s=std::sin(0.3+k);
    R(k,k)=R(k+1,k+1)=c; R(k,k+1)=-s; R(k+1,k)=s; U=U*R;
  }
  V6 lambda; lambda << 1000,500,100,10,2,0;
  const M6 original=U*lambda.asDiagonal()*U.transpose();
  const V6 residual=original*V6::Ones();
  A=original; b=residual;
  check(applySpectralReliability(A,b,o,stats), "mixed system accepted");
  check(Eigen::SelfAdjointEigenSolver<M6>(A).eigenvalues().minCoeff()>-1e-9, "PSD preserved");
  check(Eigen::SelfAdjointEigenSolver<M6>(original-A).eigenvalues().minCoeff()>-1e-9,
        "no direction gains information");
  check((b-A*V6::Ones()).norm()<1e-9, "same residual likelihood weights for Hessian and rhs");
  // Changing pose chart units while updating length scale yields the same physical system.
  M6 T=M6::Identity(); T.topLeftCorner<3,3>()*=4;
  M6 scaled=T*original*T; V6 scaled_b=T*residual;
  o.rotation_length_scale=4;
  check(applySpectralReliability(scaled,scaled_b,o,stats), "scaled chart valid");
  check((scaled-T*A*T).norm()<1e-8 && (scaled_b-T*b).norm()<1e-8, "coordinate scale covariance");
  o.rotation_length_scale=1;
  // Increasing weak-mode support to the reliable threshold releases attenuation.
  A=M6::Identity()*1000; A(5,5)=20; b=V6::Ones();
  check(applySpectralReliability(A,b,o,stats) && stats.spectral_attenuated==0,
        "recovered final geometry releases attenuation");
  A=M6::Identity(); A(0,0)=-1; b=V6::Ones();
  const M6 invalid=A;
  check(!applySpectralReliability(A,b,o,stats) && A==invalid && b==V6::Ones(), "atomic indefinite rejection");
  A.setZero(); check(!applySpectralReliability(A,b,o,stats), "empty measurement rejected");
  A.setIdentity(); b[0]=std::numeric_limits<double>::quiet_NaN();
  check(!applySpectralReliability(A,b,o,stats), "nonfinite residual rejected");
  o.enable_spectral_reliability=true;
  check(!o.valid(0.5), "spectral requires diagnostics enabled");
  o.enable_degeneracy=true; check(o.valid(0.5), "valid spectral configuration");
  o.spectral_full_ratio=0; check(!o.valid(0.5), "zero threshold rejected");
  std::cout << "spectral mathematical invariants passed\n";
}
