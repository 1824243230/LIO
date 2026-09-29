#include "lio/ESKF.h"
#include "lio/geometry.h"
#include <iostream>
#include <stdexcept>
using namespace LI2Sup;
using namespace BASIC;
void check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
int main() {
  ESKF baseline, guarded;
  auto plane = [](const ESKF::KFState& state, M6& H, V6& b) {
    H.setZero(); b.setZero();
    H(4,4)=H(5,5)=1000;
    b[4]=-1000*state.pose.t_[1]; b[5]=-1000*state.pose.t_[2];
  };
  baseline.UpdateObserve(plane);
  guarded.UpdateObserve(plane, [](const ESKF::STATE&) { return true; });
  check((baseline.GetCov()-guarded.GetCov()).norm()==0, "disabled guard preserves covariance");
  check((baseline.GetSE3().t_-guarded.GetSE3().t_).norm()==0, "disabled guard preserves state");
  // Detection-only experiment: identical nonzero plane observations over successive
  // frames, with analysis injected into only one callback. No bump rows or guard.
  ESKF plain_sequence, detected_sequence;
  geometry::Analyzer detector;
  geometry::Options detection_options;
  geometry::DegeneracyResult result;
  for(int frame=0;frame<8;++frame) {
    const scalar target=scalar(0.02*(frame+1));
    auto measurement = [&](const ESKF::KFState& state, M6& H, V6& b) {
      H.setZero(); b.setZero();
      for(int axis=0;axis<3;++axis) {
        V6 J=V6::Zero(); J[axis+3]=1;
        H += 1000*J*J.transpose();
        b -= 1000*J*(state.pose.t_[axis]-target*(axis+1));
      }
    };
    int calls=0;
    plain_sequence.UpdateObserve(measurement);
    detected_sequence.UpdateObserve([&](const ESKF::KFState& state, M6& H, V6& b) {
      measurement(state,H,b);
      const M6 before_H=H; const V6 before_b=b;
      result=detector.analyze(H.cast<double>(),detection_options,calls++==0);
      check((H-before_H).norm()==0 && (b-before_b).norm()==0, "analysis leaves measurement unchanged");
    });
    check(result.valid && result.weak_dim==3 && result.status==geometry::Status::DEGENERATE,
          "detector executes on pose information");
    const auto a=plain_sequence.GetSysState(), b=detected_sequence.GetSysState();
    check((a.R.R_-b.R.R_).norm()==0 && (a.p-b.p).norm()==0 && (a.v-b.v).norm()==0 &&
          (a.bg-b.bg).norm()==0 && (a.ba-b.ba).norm()==0 &&
          (plain_sequence.GetGravity()-detected_sequence.GetGravity()).norm()==0 &&
          (plain_sequence.GetCov()-detected_sequence.GetCov()).norm()==0,
          "detection-only preserves full filter result across frames");
  }
  check(plain_sequence.GetSE3().t_.norm()>0.01, "equivalence test has a nonzero correction");
  auto bump = [&](double variance) {
    return [&,variance](const ESKF::KFState& state, M6& H, V6& b) {
      plane(state,H,b);
      H(3,3)=1/variance;
      b[3]=(0.1-state.pose.t_[0])/variance;
    };
  };
  ESKF strong, weak;
  geometry::Options covariance_options;
  double strong_variance=0, weak_variance=0;
  check(geometry::adaptiveBumpVariance(0.9,covariance_options,strong_variance) &&
        geometry::adaptiveBumpVariance(0.05,covariance_options,weak_variance),
        "quality produces valid per-measurement covariance");
  strong.UpdateObserve(bump(strong_variance)); weak.UpdateObserve(bump(weak_variance));
  check(strong.GetSE3().t_[0]>weak.GetSE3().t_[0], "adaptive covariance changes weak-direction correction");
  check(strong.GetCov()(3,3)<weak.GetCov()(3,3), "adaptive covariance enters posterior");
  check(baseline.GetSE3().t_[0]==0, "plane leaves weak translation unobserved");
  ESKF retry, original=retry;
  bool success = retry.UpdateObserve(bump(0.001), [](const ESKF::STATE& dx) { return dx.segment<3>(3).norm()<0.001; });
  check(!success, "abnormal correction guard");
  retry=original; retry.UpdateObserve(plane);
  check((retry.GetCov()-baseline.GetCov()).norm()==0 &&
        (retry.GetSE3().t_-baseline.GetSE3().t_).norm()==0, "full snapshot restore and plane retry");
  // End-to-end equivalence to an explicit H/r/diagonal-R stack in each geometry state.
  for(auto status : {geometry::Status::NORMAL,geometry::Status::WEAK,geometry::Status::DEGENERATE}) {
    for(bool useful_bumps : {false,true}) {
      ESKF information_filter, explicit_filter;
      geometry::Options opts;
      geometry::DegeneracyResult mode;
      mode.valid=true; mode.status=status; mode.ratio=status==geometry::Status::DEGENERATE ? 0 : 0.009;
      mode.weak[3]=true; mode.weak_dim=1; // Current weak translation direction is world x.
      auto observe_stack = [&](bool explicit_stack,const ESKF::KFState& state,M6& A,V6& b) {
        std::vector<geometry::V6> Jp(4,geometry::V6::Zero());
        Jp[0][4]=1; Jp[1][5]=1; Jp[2][5]=1; Jp[3][3]=1e6;
        const std::vector<unsigned char> valid={1,1,1,0}; // Fourth plane failed reliability gate.
        const std::vector<double> rp={double(state.pose.t_[1])-0.02,
          double(state.pose.t_[2])-0.01,double(state.pose.t_[2])+0.02,100};
        geometry::M6 plane_A=geometry::M6::Zero(); geometry::V6 plane_b=geometry::V6::Zero();
        for(size_t i=0;i<Jp.size();++i) if(valid[i]) {
          plane_A+=1000*Jp[i]*Jp[i].transpose(); plane_b-=1000*Jp[i]*rp[i];
        }
        std::vector<geometry::Candidate> selected;
        geometry::Diagnostics stats;
        for(size_t i : {size_t(1),size_t(2)}) {
          geometry::Surface surface;
          surface.normal=geometry::V3(i==1 ? 0.4 : -0.3,0,1);
          surface.residual=surface.normal.dot(state.pose.t_.cast<double>())+(i==1 ? -0.01 : 0.02);
          surface.mid=opts.bump_mid_scale;
          surface.gradient=useful_bumps ? std::abs(surface.normal.x()) : 0;
          surface.confidence=opts.bump_pixel_scale*(i==1 ? 1.0 : 0.5);
          geometry::Candidate c; c.index=i;
          if(geometry::candidate(surface,geometry::V3::Zero(),geometry::M3::Identity(),mode,opts,c,stats)) selected.push_back(c);
        }
        geometry::selectWeakDirectionConstraints(selected,mode,opts);
        check(selected.size()==size_t(useful_bumps && status!=geometry::Status::NORMAL ? 2 : 0),
              "state-specific stack activates only reliable degenerate bump rows");
        // Build the independent reference from only retained plane rows plus selected bump rows.
        Eigen::Matrix<double,Eigen::Dynamic,6> H(3,6);
        Eigen::VectorXd residual(3), variance(3);
        int row=0;
        for(size_t i=0;i<Jp.size();++i) if(valid[i]) {
          bool replaced=false;
          for(const auto& c:selected) if(c.index==i) replaced=true;
          if(replaced) continue;
          H.row(row)=Jp[i].transpose(); residual[row]=rp[i]; variance[row]=0.001; ++row;
        }
        for(const auto& c:selected) {
          H.row(row)=c.J.transpose(); residual[row]=c.residual; variance[row]=c.variance; ++row;
        }
        check(row==3,"mutually exclusive stack preserves one row per accepted point");
        const geometry::M6 reference_A=H.transpose()*variance.cwiseInverse().asDiagonal()*H;
        const geometry::V6 reference_b=-H.transpose()*variance.cwiseInverse().asDiagonal()*residual;
        geometry::applyExclusiveBumps(Jp,rp,valid,selected,plane_A,plane_b);
        check((plane_A-reference_A).norm()<1e-8 && (plane_b-reference_b).norm()<1e-10,
              "information update equals explicit stacked H/r/R");
        check(std::abs(plane_A(4,4)-1000)<1e-10,"degeneracy never globally increases retained plane precision");
        A=(explicit_stack ? reference_A : plane_A).cast<scalar>();
        b=(explicit_stack ? reference_b : plane_b).cast<scalar>();
      };
      information_filter.UpdateObserve([&](const ESKF::KFState& state,M6& A,V6& b) { observe_stack(false,state,A,b); });
      explicit_filter.UpdateObserve([&](const ESKF::KFState& state,M6& A,V6& b) { observe_stack(true,state,A,b); });
      const auto a=information_filter.GetSysState(), b=explicit_filter.GetSysState();
      check((a.R.R_-b.R.R_).norm()<1e-6 && (a.p-b.p).norm()<1e-6 && (a.v-b.v).norm()<1e-6 &&
            (a.bg-b.bg).norm()<1e-6 && (a.ba-b.ba).norm()<1e-6 &&
            (information_filter.GetGravity()-explicit_filter.GetGravity()).norm()<1e-6 &&
            (information_filter.GetCov()-explicit_filter.GetCov()).norm()<1e-6,
            "all-state posterior matches explicit measurement stack");
    }
  }
  std::cout << "filter tests passed: baseline equivalence, adaptive information, guarded rollback\n";
}
