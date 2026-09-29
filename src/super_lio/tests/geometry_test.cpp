#include "lio/geometry.h"
#include <Eigen/Geometry>
#include <iostream>
#include <stdexcept>
#include <limits>
using namespace LI2Sup::geometry;
void check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
int main() {
  Options o;
  check(o.valid(0.5), "default options");
  Options bad = o; bad.sigma_b_min = -1; check(!bad.valid(0.5), "invalid variance");
  Analyzer analyzer;
  M6 G = M6::Identity();
  check(analyzer.analyze(G,o,true).status == Status::NORMAL, "normal geometry");
  G(5,5) = 0.005;
  check(analyzer.analyze(G,o,true).status == Status::WEAK, "enter weak");
  G(5,5) = 0.015;
  check(analyzer.analyze(G,o,true).status == Status::WEAK, "weak hysteresis");
  G(5,5) = 0.0001;
  auto d = analyzer.analyze(G,o,true);
  check(d.status == Status::DEGENERATE && d.weak_dim == 1, "degenerate eigenspace");
  G(5,5) = 0.002;
  check(analyzer.analyze(G,o,true).status == Status::DEGENERATE, "degenerate hysteresis");
  G(5,5) = 0.03;
  check(analyzer.analyze(G,o,true).status == Status::NORMAL, "hysteresis exit");
  check(!analyzer.analyze(M6::Zero(),o,true).valid, "empty information fallback");
  G(0,0) = std::numeric_limits<double>::quiet_NaN();
  check(!analyzer.analyze(G,o,true).valid, "NaN information fallback");

  G = M6::Identity(); G(0,0) = std::numeric_limits<double>::infinity();
  check(!analyzer.analyze(G,o,true).valid, "Inf information fallback");
  check(!analyzer.analyze(1e-14*M6::Identity(),o,true).valid, "near-zero maximum eigenvalue");
  check(!analyzer.analyze(std::numeric_limits<double>::max()*M6::Identity(),o,true).valid,
        "symmetrization overflow fallback");
  G = M6::Identity(); G(0,0)=-0.1;
  check(!analyzer.analyze(G,o,true).valid, "indefinite information fallback");
  Analyzer jitter;
  G = M6::Identity(); G(5,5)=0.009;
  auto ordered=jitter.analyze(G,o,true);
  check(ordered.status==Status::WEAK, "jitter starts weak");
  for(int i=0;i<20;++i) {
    G(5,5)=i%2 ? 0.009 : 0.011;
    check(jitter.analyze(G,o,true).status==Status::WEAK, "threshold jitter does not toggle state");
  }
  G(5,5)=0.1;
  check(jitter.analyze(G,o,false).status==Status::WEAK, "iteration does not advance hysteresis");
  check(jitter.analyze(G,o,true).status==Status::NORMAL, "next frame advances hysteresis");
  for(int k=0;k<5;++k) check(ordered.eigenvalues[k]>=ordered.eigenvalues[k+1], "descending spectrum");
  int weak_count=0;
  for(int k=0;k<6;++k) weak_count += ordered.eigenvalues[k]/(ordered.eigenvalues[0]+1e-12)<o.degeneracy_eigen_ratio_threshold;
  check(weak_count==ordered.weak_dim, "weak dimension matches spectrum");
  check(std::abs(ordered.ratio-ordered.eigenvalues[5]/(ordered.eigenvalues[0]+1e-12))<1e-15,
        "ratio matches regularized spectrum");

  // Nontrivial world/local chart and right body rotation; central difference all six columns.
  const M3 R = Eigen::AngleAxisd(0.6,V3(1,2,3).normalized()).toRotationMatrix();
  const M3 C = Eigen::AngleAxisd(-0.4,V3(2,-1,1).normalized()).toRotationMatrix();
  const V3 body(0.2,0.4,0.7), t(1,2,3);
  const V3 normal = C.transpose()*V3(-0.2,0.3,1);
  const V6 J = poseJacobian(body,R,normal);
  auto residual = [&](const M3& rotation, const V3& trans) {
    V3 p = C*(rotation*body+trans); return p.z()-0.2*p.x()+0.3*p.y();
  };
  for (int k=0;k<6;++k) {
    const double eps = 1e-6;
    M3 plus=R, minus=R; V3 tp=t, tm=t;
    if (k<3) {
      plus = R*Eigen::AngleAxisd(eps,V3::Unit(k)).toRotationMatrix();
      minus = R*Eigen::AngleAxisd(-eps,V3::Unit(k)).toRotationMatrix();
    } else { tp[k-3]+=eps; tm[k-3]-=eps; }
    check(std::abs((residual(plus,tp)-residual(minus,tm))/(2*eps)-J[k])<1e-8, "right perturbation Jacobian");
  }
  // A flat chart contributes no tangent translation information. A height gradient
  // changes the measurement direction, but a single scalar row still has rank one.
  const V3 chart_normal = C.transpose()*V3::UnitZ();
  const V3 tangent_x = C.transpose()*V3::UnitX();
  const V3 tangent_y = C.transpose()*V3::UnitY();
  const V6 flat = poseJacobian(body,R,chart_normal);
  check(std::abs(flat.tail<3>().dot(tangent_x))<1e-12 &&
        std::abs(flat.tail<3>().dot(tangent_y))<1e-12, "flat chart tangent nullspace");
  check(std::abs(J.tail<3>().dot(tangent_x)+0.2)<1e-12 &&
        std::abs(J.tail<3>().dot(tangent_y)-0.3)<1e-12, "height gradient tangent information");
  Options scaled_options = o; scaled_options.rotation_length_scale=3;
  M6 scaling = M6::Identity(); scaling.topLeftCorner<3,3>() /= 3;
  const M6 information = M6::Identity()+J*J.transpose();
  Analyzer scaled_analyzer;
  const auto scaled_result = scaled_analyzer.analyze(information,scaled_options,true);
  check(scaled_result.valid &&
        (scaled_result.eigenvectors*scaled_result.eigenvalues.asDiagonal()*scaled_result.eigenvectors.transpose()-
         scaling*information*scaling).norm()<1e-10, "scaled geometry eigendecomposition");
  Surface s; s.normal=V3(0,0,1); s.mid=0.02; s.gradient=0.3; s.confidence=2;
  Candidate candidate_good; Diagnostics stats;
  check(candidate(s,V3::Zero(),M3::Identity(),d,o,candidate_good,stats), "weak contribution accepted");
  check(candidate_good.variance>=o.sigma_b_min*o.sigma_b_min && candidate_good.variance<=o.sigma_b_max*o.sigma_b_max, "covariance bounds");
  s.confidence=0.5;
  Candidate candidate_low;
  check(candidate(s,V3::Zero(),M3::Identity(),d,o,candidate_low,stats), "lower confidence accepted");
  check(candidate_low.variance>candidate_good.variance, "adaptive covariance trust");
  s.normal=V3(1,0,0);
  check(!candidate(s,V3::Zero(),M3::Identity(),d,o,candidate_low,stats), "unhelpful direction rejected");
  s.normal=V3(0,0,1); s.gradient=0;
  check(!candidate(s,V3::Zero(),M3::Identity(),d,o,candidate_low,stats), "flat surface rejected");
  s.gradient=0.3; s.residual=1;
  check(!candidate(s,V3::Zero(),M3::Identity(),d,o,candidate_low,stats), "outlier rejected");

  Surface half;
  half.mid=o.bump_mid_scale*0.5; half.gradient=o.bump_gradient_scale*0.5;
  half.confidence=o.bump_pixel_scale*0.5;
  half.normal=V3(0,0,std::sqrt(o.bump_weak_scale*0.5));
  Candidate quality_test;
  check(candidate(half,V3::Zero(),M3::Identity(),d,o,quality_test,stats), "fractional quality accepted");
  check(std::abs(quality_test.w_mid-0.5)<1e-12 && std::abs(quality_test.w_grad-0.5)<1e-12 &&
        std::abs(quality_test.w_pixel-0.5)<1e-12 && std::abs(quality_test.w_weak-0.5)<1e-12 &&
        std::abs(quality_test.local_quality-0.0625)<1e-12 && std::abs(quality_test.quality-0.0625)<1e-12,
        "four clamped quality factors multiply exactly");
  half.residual=2*o.bump_huber_delta;
  check(candidate(half,V3::Zero(),M3::Identity(),d,o,quality_test,stats) &&
        std::abs(quality_test.robust_weight-0.5)<1e-12 && std::abs(quality_test.quality-0.03125)<1e-12,
        "Huber is separate from local quality and cannot amplify it");
  half.residual=0; half.mid=half.gradient=half.confidence=std::numeric_limits<double>::max();
  half.normal=V3(0,0,1);
  check(candidate(half,V3::Zero(),M3::Identity(),d,o,quality_test,stats) &&
        quality_test.w_mid==1 && quality_test.w_grad==1 && quality_test.w_pixel==1 &&
        quality_test.w_weak==1 && quality_test.quality==1, "all quality factors saturate at one");
  half.normal=V3(0,0,1e200);
  check(!candidate(half,V3::Zero(),M3::Identity(),d,o,quality_test,stats), "overflowed weak energy cannot become high confidence");
  half.normal=V3(0,0,1); half.mid=-1;
  check(!candidate(half,V3::Zero(),M3::Identity(),d,o,quality_test,stats), "negative MID rejected");
  half.mid=std::numeric_limits<double>::quiet_NaN();
  check(!candidate(half,V3::Zero(),M3::Identity(),d,o,quality_test,stats), "NaN quality input rejected");

  // MID is a per-observed-pixel mean, not a confidence-weighted mean or a full-image mean.
  BumpLayer mid_layer;
  mid_layer.width=4; mid_layer.height=2;
  mid_layer.image={-0.02,0.04,100.0,0.0,std::numeric_limits<double>::quiet_NaN(),
                   std::numeric_limits<double>::infinity(),0.8,0.9};
  mid_layer.weight={1,100,0,2,1,1,std::numeric_limits<double>::infinity(),-1};
  mid_layer.updateMid();
  check(mid_layer.valid && std::abs(mid_layer.mid-0.02)<1e-14,
        "MID averages absolute finite observed pixels including observed zero");
  mid_layer.weight[1]=100000;
  mid_layer.updateMid();
  check(std::abs(mid_layer.mid-0.02)<1e-14, "pixel confidence does not reweight MID");
  std::fill(mid_layer.weight.begin(),mid_layer.weight.end(),0);
  mid_layer.updateMid();
  check(!mid_layer.valid && mid_layer.mid==0, "empty observed set invalidates layer");
  mid_layer.weight[0]=1; mid_layer.image[0]=0;
  mid_layer.updateMid();
  check(mid_layer.valid && mid_layer.mid==0, "observed flat zero is a valid pixel");
  mid_layer.weight.pop_back(); mid_layer.updateMid();
  check(!mid_layer.valid, "malformed image dimensions rejected");

  s.mid=100; s.residual=0; s.gradient=0; s.confidence=2; s.normal=V3(0,0,1);
  check(!candidate(s,V3::Zero(),M3::Identity(),d,o,candidate_low,stats), "high MID cannot bypass gradient gate");
  s.gradient=0.3; s.confidence=0;
  check(!candidate(s,V3::Zero(),M3::Identity(),d,o,candidate_low,stats), "high MID cannot bypass pixel gate");
  s.confidence=2; s.normal=V3(1,0,0);
  check(!candidate(s,V3::Zero(),M3::Identity(),d,o,candidate_low,stats), "high MID cannot bypass weak-direction gate");

  // Mutually exclusive fusion must equal a directly assembled mixed measurement stack.
  std::vector<V6> plane_rows(3,V6::Zero());
  plane_rows[0][3]=1; plane_rows[1][4]=1; plane_rows[2][5]=1;
  const std::vector<double> residuals={0.1,-0.2,0.3};
  const std::vector<unsigned char> valid_planes={1,1,1};
  M6 original_info=M6::Zero(); V6 original_rhs=V6::Zero();
  for(size_t i=0;i<plane_rows.size();++i) {
    original_info+=1000*plane_rows[i]*plane_rows[i].transpose();
    original_rhs-=1000*plane_rows[i]*residuals[i];
  }
  Candidate replacement; replacement.index=1; replacement.J=plane_rows[1]+0.5*plane_rows[0];
  replacement.residual=0.05; replacement.variance=0.002;
  M6 mixed_info=original_info; V6 mixed_rhs=original_rhs;
  std::vector<Candidate> selected={replacement,replacement};
  applyExclusiveBumps(plane_rows,residuals,valid_planes,selected,mixed_info,mixed_rhs);
  M6 expected_info=M6::Zero(); V6 expected_rhs=V6::Zero();
  for(size_t i : {size_t(0),size_t(2)}) {
    expected_info+=1000*plane_rows[i]*plane_rows[i].transpose();
    expected_rhs-=1000*plane_rows[i]*residuals[i];
  }
  expected_info+=replacement.J*replacement.J.transpose()/replacement.variance;
  expected_rhs-=replacement.J*replacement.residual/replacement.variance;
  check(selected.size()==1 && (mixed_info-expected_info).norm()<1e-12 &&
        (mixed_rhs-expected_rhs).norm()<1e-12,"exclusive stack removes exactly one plane row and deduplicates");
  // Identical bump/plane rows at identical precision must not double their information.
  replacement.J=plane_rows[1]; replacement.residual=residuals[1]; replacement.variance=0.001;
  mixed_info=original_info; mixed_rhs=original_rhs; selected={replacement};
  applyExclusiveBumps(plane_rows,residuals,valid_planes,selected,mixed_info,mixed_rhs);
  check((mixed_info-original_info).norm()<1e-12 && (mixed_rhs-original_rhs).norm()<1e-12,
        "identical replacement does not double count plane information");
  replacement.variance=-1; selected={replacement};
  applyExclusiveBumps(plane_rows,residuals,valid_planes,selected,mixed_info,mixed_rhs);
  check(selected.empty() && (mixed_info-original_info).norm()<1e-12 && (mixed_rhs-original_rhs).norm()<1e-12,
        "rejected covariance keeps original plane");
  replacement.variance=0.001; replacement.index=100; selected={replacement};
  applyExclusiveBumps(plane_rows,residuals,valid_planes,selected,mixed_info,mixed_rhs);
  check(selected.empty() && (mixed_info-original_info).norm()==0,"unknown point cannot remove a plane row");
  applyExclusiveBumps(plane_rows,residuals,valid_planes,selected,mixed_info,mixed_rhs);
  check((mixed_info-original_info).norm()==0 && (mixed_rhs-original_rhs).norm()==0,
        "no selected bump preserves pure plane measurements");

  // Invalid bump rows must neither erase their plane nor block subsequent valid rows.
  for (int failure=0;failure<5;++failure) {
    Candidate bad_row=replacement;
    bad_row.index=0; bad_row.J=plane_rows[0]; bad_row.variance=0.001;
    if(failure==0) bad_row.J[0]=std::numeric_limits<double>::quiet_NaN();
    if(failure==1) bad_row.J[0]=std::numeric_limits<double>::infinity();
    if(failure==2) bad_row.variance=std::numeric_limits<double>::quiet_NaN();
    if(failure==3) bad_row.variance=std::numeric_limits<double>::infinity();
    if(failure==4) bad_row.variance=0;
    Candidate good_row=bad_row;
    good_row.index=1; good_row.J=plane_rows[1]; good_row.residual=residuals[1]; good_row.variance=0.001;
    selected={bad_row,good_row}; mixed_info=original_info; mixed_rhs=original_rhs;
    applyExclusiveBumps(plane_rows,residuals,valid_planes,selected,mixed_info,mixed_rhs);
    check(selected.size()==1 && selected[0].index==1 &&
          (mixed_info-original_info).norm()==0 && (mixed_rhs-original_rhs).norm()==0,
          "invalid Jacobian/covariance preserves plane and allows subsequent valid measurement");
  }

  // Oblique weak eigendirection: large MID alone must not outrank directional usefulness.
  DegeneracyResult oblique; oblique.valid=true; oblique.status=Status::WEAK;
  oblique.weak_dim=1; oblique.weak[0]=true;
  oblique.eigenvectors.col(0) << 0,0,0,1/std::sqrt(2.0),1/std::sqrt(2.0),0;
  Surface strong_only; strong_only.mid=1; strong_only.gradient=1; strong_only.confidence=10;
  strong_only.normal=V3(1,-1,0).normalized();
  Candidate direction_test;
  check(!candidate(strong_only,V3::Zero(),M3::Identity(),oblique,o,direction_test,stats),
        "high-MID observation orthogonal to oblique weak mode rejected");
  strong_only.mid=o.bump_mid_scale*0.5; strong_only.normal=V3(1,1,0).normalized();
  check(candidate(strong_only,V3::Zero(),M3::Identity(),oblique,o,direction_test,stats) &&
        std::abs(direction_test.weak_score-1)<1e-12,
        "lower-MID observation aligned with oblique weak mode accepted");
  Options select_options=o; select_options.max_bump_constraints=2; select_options.weak_direction_score_threshold=0.2;
  auto ranked_candidate=[&](size_t index,double score) {
    Candidate c=direction_test; c.index=index; c.quality=score; return c;
  };
  std::vector<Candidate> ranked={ranked_candidate(5,0.7),ranked_candidate(2,0.8),ranked_candidate(2,0.9),
                                ranked_candidate(1,0.7),ranked_candidate(8,0.2),ranked_candidate(9,0.1)};
  const auto original_ranked=ranked;
  selectWeakDirectionConstraints(ranked,oblique,select_options);
  check(ranked.size()==2 && ranked[0].index==2 && ranked[0].quality==0.9 && ranked[1].index==1,
        "Top-K uses quality, distinct indices, deterministic ties");
  ranked=original_ranked; select_options.max_bump_constraints=20;
  selectWeakDirectionConstraints(ranked,oblique,select_options);
  check(ranked.size()==3 && ranked.back().index==5,"strict threshold excludes equality, cap does not add observations");
  auto normal_selection=oblique; normal_selection.status=Status::NORMAL;
  selectWeakDirectionConstraints(ranked,normal_selection,select_options);
  check(ranked.empty(),"NORMAL clears final bump selection");
  ranked=original_ranked;
  selectWeakDirectionConstraints(ranked,DegeneracyResult{},select_options);
  check(ranked.empty(),"invalid geometry clears final bump selection");

  double previous_variance=std::numeric_limits<double>::infinity(), adaptive_variance=0;
  for(double q : {0.0,0.001,0.01,0.1,0.5,1.0}) {
    check(adaptiveBumpVariance(q,o,adaptive_variance), "valid adaptive covariance");
    const double expected=std::clamp(o.sigma_b_base*o.sigma_b_base/(q+1e-12),
                                    o.sigma_b_min*o.sigma_b_min,o.sigma_b_max*o.sigma_b_max);
    check(std::abs(adaptive_variance-expected)<1e-14 && adaptive_variance<=previous_variance,
          "base covariance formula and monotonic trust");
    previous_variance=adaptive_variance;
  }
  check(adaptiveBumpVariance(0,o,adaptive_variance) && adaptive_variance==o.sigma_b_max*o.sigma_b_max,
        "low quality saturates at maximum variance");
  check(adaptiveBumpVariance(1,o,adaptive_variance,1e-4) && adaptive_variance==o.sigma_b_min*o.sigma_b_min,
        "severity cannot cross minimum variance");
  double base_variance=0, severity_variance=0;
  check(adaptiveBumpVariance(0.5,o,base_variance) && adaptiveBumpVariance(0.5,o,severity_variance,0.5) &&
        std::abs(severity_variance/base_variance-0.5)<1e-12, "optional severity factor scales numerator");
  check(!adaptiveBumpVariance(-0.1,o,adaptive_variance) &&
        !adaptiveBumpVariance(1.1,o,adaptive_variance) &&
        !adaptiveBumpVariance(std::numeric_limits<double>::quiet_NaN(),o,adaptive_variance),
        "invalid quality rejects covariance");
  Options invalid_sigma=o; invalid_sigma.sigma_b_base=std::numeric_limits<double>::infinity();
  check(!adaptiveBumpVariance(0.5,invalid_sigma,adaptive_variance), "nonfinite sigma rejected");
  invalid_sigma=o; invalid_sigma.sigma_b_min=0;
  check(!adaptiveBumpVariance(0.5,invalid_sigma,adaptive_variance), "nonpositive variance bound rejected");

  double severity=0;
  check(degeneracySeverity(o.degeneracy_ratio_exit_weak,o,severity) && severity==0,"normal severity endpoint");
  check(degeneracySeverity(0,o,severity) && severity==1,"strong degeneracy severity endpoint");
  check(degeneracySeverity((o.degeneracy_ratio_exit_weak+o.degeneracy_ratio_enter_degenerate)/2,o,severity) &&
        std::abs(severity-0.5)<1e-9,"severity interpolation midpoint");
  check(degeneracySeverity(1,o,severity) && severity==0,"severity lower clamp");
  check(!degeneracySeverity(-0.1,o,severity) &&
        !degeneracySeverity(std::numeric_limits<double>::quiet_NaN(),o,severity),"invalid ratio rejected");
  Options bad_severity=o; bad_severity.degeneracy_ratio_exit_weak=o.degeneracy_ratio_enter_degenerate;
  check(!degeneracySeverity(0,bad_severity,severity),"invalid severity interval rejected");
  Surface reliable; reliable.mid=o.bump_mid_scale; reliable.gradient=o.bump_gradient_scale;
  reliable.confidence=o.bump_pixel_scale; reliable.normal=V3(0,0,1);
  auto global_state=d;
  Candidate global_candidate;
  double prior_R=std::numeric_limits<double>::infinity();
  for(double ratio : {0.015,0.009,0.003,0.001,0.0}) {
    global_state.ratio=ratio;
    global_state.status=ratio>=0.001 ? Status::WEAK : Status::DEGENERATE;
    check(candidate(reliable,V3::Zero(),M3::Identity(),global_state,o,global_candidate,stats),
          "reliable candidate survives severity sweep");
    check(global_candidate.local_quality==1 && global_candidate.quality==1 &&
          global_candidate.degeneracy_severity>=0 && global_candidate.degeneracy_severity<=1 &&
          global_candidate.variance<=prior_R,"severity reduces variance without changing local quality");
    prior_R=global_candidate.variance;
  }
  check(std::abs(global_candidate.variance_factor-o.bump_degenerate_variance_factor)<1e-12,
        "severe variance approaches configured lower base variance");
  for(int gate=0;gate<4;++gate) {
    Surface unreliable=reliable;
    if(gate==0) unreliable.mid=0;
    if(gate==1) unreliable.gradient=0;
    if(gate==2) unreliable.confidence=0;
    if(gate==3) unreliable.normal=V3(1,0,0);
    check(!candidate(unreliable,V3::Zero(),M3::Identity(),global_state,o,global_candidate,stats),
          "maximal severity cannot bypass any local reliability gate");
  }
  Options no_severity=o; no_severity.bump_degenerate_variance_factor=1;
  global_state.ratio=0.015;
  check(candidate(reliable,V3::Zero(),M3::Identity(),global_state,no_severity,global_candidate,stats),"severity ablation weak");
  const double unchanged_R=global_candidate.variance;
  global_state.ratio=0;
  check(candidate(reliable,V3::Zero(),M3::Identity(),global_state,no_severity,global_candidate,stats) &&
        global_candidate.variance==unchanged_R,"factor one disables global severity modulation");

  // Bilinear safety gate: each of four real pixels is mandatory, including image borders.
  BumpLayer tile; tile.width=2; tile.height=2; tile.resolution=0.1;
  tile.image={0.01,0.02,0.03,0.04}; tile.weight={1,2,3,4}; tile.updateMid();
  Surface tile_result;
  check(tile.query(V3(0,0,0.025),tile_result) && std::abs(tile_result.residual)<1e-12 &&
        tile_result.confidence==1,"valid four-pixel interpolation");
  for(int i=0;i<4;++i) {
    auto missing=tile; missing.weight[i]=0;
    check(!missing.query(V3::Zero(),tile_result) && tile_result.confidence==0,"each missing bilinear pixel rejects and clears output");
  }
  check(!tile.query(V3(0.05,0,0),tile_result),"last pixel column lacks full bilinear support");
  check(!tile.query(V3(0,0.05,0),tile_result),"last pixel row lacks full bilinear support");
  auto broken=tile; broken.valid=false;
  check(!broken.query(V3::Zero(),tile_result),"invalid layer rejected");
  broken=tile; broken.image.pop_back();
  check(!broken.query(V3::Zero(),tile_result),"inconsistent layer arrays rejected");
  for(double nonfinite : {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
    broken=tile; broken.image[0]=nonfinite;
    check(!broken.query(V3::Zero(),tile_result),"nonfinite observed height rejected");
    broken=tile; broken.weight[0]=nonfinite;
    check(!broken.query(V3::Zero(),tile_result),"nonfinite pixel confidence rejected");
    broken=tile; broken.R_CG(0,0)=nonfinite;
    check(!broken.query(V3::Zero(),tile_result),"nonfinite chart rejected");
    check(!tile.query(V3(nonfinite,0,0),tile_result),"nonfinite query point rejected");
    for(int field=0;field<5;++field) {
      Surface invalid=reliable;
      if(field==0) invalid.mid=nonfinite;
      if(field==1) invalid.gradient=nonfinite;
      if(field==2) invalid.confidence=nonfinite;
      if(field==3) invalid.residual=nonfinite;
      if(field==4) invalid.normal[0]=nonfinite;
      check(!candidate(invalid,V3::Zero(),M3::Identity(),d,o,global_candidate,stats),"nonfinite measurement field rejected");
    }
  }
  broken=tile; broken.resolution=std::numeric_limits<double>::denorm_min();
  check(!broken.query(V3(1,0,0),tile_result),"coordinate overflow rejected before integer conversion");
  reliable.residual=std::nextafter(o.bump_residual_max,0.0);
  check(candidate(reliable,V3::Zero(),M3::Identity(),d,o,global_candidate,stats),"residual strictly inside boundary accepted");
  reliable.residual=-reliable.residual;
  check(candidate(reliable,V3::Zero(),M3::Identity(),d,o,global_candidate,stats),"negative residual strictly inside boundary accepted");
  reliable.residual=o.bump_residual_max;
  check(!candidate(reliable,V3::Zero(),M3::Identity(),d,o,global_candidate,stats),"strict positive residual boundary rejected");
  reliable.residual=-o.bump_residual_max;
  check(!candidate(reliable,V3::Zero(),M3::Identity(),d,o,global_candidate,stats),"strict negative residual boundary rejected");
  reliable.residual=std::nextafter(o.bump_residual_max,std::numeric_limits<double>::infinity());
  check(!candidate(reliable,V3::Zero(),M3::Identity(),d,o,global_candidate,stats),"positive residual above limit rejected");
  reliable.residual=-reliable.residual;
  check(!candidate(reliable,V3::Zero(),M3::Identity(),d,o,global_candidate,stats),"negative residual above absolute limit rejected");
  check(global_candidate.variance==0 && global_candidate.quality==0,"rejected candidate clears prior acceptance data");

  const double delta=o.bump_huber_delta;
  for(double r : {0.0,delta/2,delta,-delta/2,-delta})
    check(huberWeight(r,delta)==1,"Huber unit weight inside and at threshold");
  check(std::abs(huberWeight(2*delta,delta)-0.5)<1e-15 &&
        huberWeight(-2*delta,delta)==huberWeight(2*delta,delta),"Huber tail weight and sign symmetry");
  check(huberWeight(0,1e-15)==1 && huberWeight(1e-15,1e-15)==1,"small valid delta retains unit inlier weight");
  check(huberWeight(1,0)==0 && huberWeight(1,-1)==0 &&
        huberWeight(std::numeric_limits<double>::infinity(),delta)==0 &&
        huberWeight(0,std::numeric_limits<double>::quiet_NaN())==0,"invalid Huber inputs reject");
  Surface kernel_surface; kernel_surface.mid=o.bump_mid_scale; kernel_surface.gradient=o.bump_gradient_scale;
  kernel_surface.confidence=o.bump_pixel_scale; kernel_surface.normal=V3(0,0,1);
  Options kernel_options=o; kernel_options.bump_degenerate_variance_factor=1;
  for(bool adaptive : {false,true}) {
    kernel_options.enable_adaptive_covariance=adaptive;
    Candidate inlier,outlier,negative;
    kernel_surface.residual=delta;
    check(candidate(kernel_surface,V3::Zero(),M3::Identity(),d,kernel_options,inlier,stats),"Huber boundary candidate");
    kernel_surface.residual=2*delta;
    check(candidate(kernel_surface,V3::Zero(),M3::Identity(),d,kernel_options,outlier,stats),"Huber downweighted candidate");
    kernel_surface.residual=-2*delta;
    check(candidate(kernel_surface,V3::Zero(),M3::Identity(),d,kernel_options,negative,stats),"negative Huber candidate");
    check(inlier.local_quality==outlier.local_quality && outlier.quality==0.5*inlier.quality &&
          std::abs(outlier.variance/inlier.variance-2)<1e-9 && negative.variance==outlier.variance,
          "Huber enters quality and covariance once in adaptive and fixed modes");
  }
  kernel_options.enable_adaptive_covariance=true;
  kernel_options.weak_direction_score_threshold=0.75;
  check(!candidate(kernel_surface,V3::Zero(),M3::Identity(),d,kernel_options,global_candidate,stats) &&
        global_candidate.variance==0,"Huber-degraded quality may reject before covariance");

  // Dense synthetic ripple, enough independent samples to populate bilinear neighbors.
  o.bump_history_points=4096; o.bump_min_roughness=0.0001; o.bump_resolution=0.04;
  BumpMap map(o,0.5);
  std::vector<V3> points;
  for(int x=0;x<40;++x) for(int y=0;y<40;++y) {
    double px=0.05+x*0.01, py=0.05+y*0.01;
    points.emplace_back(px,py,0.2+0.015*std::sin(40*px)*std::cos(30*py));
  }
  map.insert(points,V3::Zero());
  Surface query;
  check(!map.query(V3(0.22,0.23,0.2),query), "unstable new chart gated");
  map.insert(points,V3::Zero());
  V3 p(0.2213,0.2337,0.2);
  check(map.query(p,query) && query.mid>0 && query.gradient>0, "ripple query");
  for(int k=0;k<3;++k) {
    Surface a,b; V3 pa=p,pb=p; pa[k]+=1e-6; pb[k]-=1e-6;
    check(map.query(pa,a) && map.query(pb,b), "bilinear neighbors");
    check(std::abs((a.residual-b.residual)/2e-6-query.normal[k])<1e-6,"bilinear gradient derivative");
  }
  // End-to-end chain: right pose perturbation -> actual bilinear map residual.
  const V3 mapped_body = R.transpose()*(p-t);
  const V6 map_J = poseJacobian(mapped_body,R,query.normal);
  for(int k=0;k<6;++k) {
    const double eps=1e-7;
    M3 rp=R, rm=R; V3 tp=t, tm=t;
    if(k<3) {
      rp=R*Eigen::AngleAxisd(eps,V3::Unit(k)).toRotationMatrix();
      rm=R*Eigen::AngleAxisd(-eps,V3::Unit(k)).toRotationMatrix();
    } else { tp[k-3]+=eps; tm[k-3]-=eps; }
    Surface a,b;
    check(map.query(rp*mapped_body+tp,a) && map.query(rm*mapped_body+tm,b), "pose query neighbors");
    check(std::abs((a.residual-b.residual)/(2*eps)-map_J[k])<1e-6, "full bump pose Jacobian");
  }
  check(!map.query(V3(5,5,5),query), "unobserved query fallback");
  Options bounded=o; bounded.bump_capacity=2;
  BumpMap small(bounded,0.5);
  small.insert({V3(0,0,0),V3(1,0,0),V3(2,0,0)},V3::Zero());
  check(small.size()==2, "bounded map eviction");
  std::cout << "geometry tests passed: hysteresis, invalid input, Jacobian, selection, covariance, height map, eviction\n";
}
