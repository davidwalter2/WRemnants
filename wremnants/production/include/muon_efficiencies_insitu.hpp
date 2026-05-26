#ifndef WREMNANTS_MUON_EFFICIENCIES_INSITU_H
#define WREMNANTS_MUON_EFFICIENCIES_INSITU_H

#include <array>
#include <boost/histogram/axis.hpp>
#include <cmath>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <memory>

namespace wrem {

// In-situ muon efficiency helper for the Z->mumu dilepton analysis.
//
// Floats the ID (idip), HLT (trigger) and Iso efficiency data/MC scale
// factors as (NCoeff-1)-order Chebyshev polynomials in muon pt, decorrelated
// per NEta-bin probe-eta (and per charge for idip). The scale factor of step
// X for a leg is
//   SF_X = exp( sum_{k=0..NCoeff-1} theta_{X,etaBin,[q],k} * T_k(xtilde) )
// with xtilde = 2*(clamp(pt,ptmin,ptmax)-ptmin)/(ptmax-ptmin) - 1.
// A leg PASSING step X contributes a factor SF_X; a leg FAILING contributes
//   (1 - SF_X * eMC) / (1 - eMC),   eMC = effMC_X(pt, eta, q).
// At theta = 0 every factor equals 1, so the nominal event weight is
// unchanged automatically.
//
// rabbit treats every output-tensor bin as an independent linearised
// (log-normal) nuisance around the nominal. We therefore store, per
// coefficient c, the analytic first-order response of the event weight:
//   res(c) = w_nom * exp( delta * d lnW / d theta_c )
// with d ln f / d s = +1 (pass) or -eMC/(1-eMC) (fail). The coefficients are
// unconstrained, so "delta" only sets the units of the nuisance: the fitted
// nuisance n_c corresponds to a Chebyshev coefficient theta_c = delta * n_c
// (physical results are delta independent). delta is kept small because the
// fail-leg derivative -eMC/(1-eMC) is large for high-efficiency steps
// (eMC ~ 0.99 for idip), which would otherwise overflow exp(). effMC is also
// clamped to [0, effMC_max] to keep the derivative finite.
// Bins that no leg/step touches stay at w_nom (exp(0) = 1).
//
// category: 0 = nominal, 1 = failIso, 2 = failHLT, 3 = failID. The tag leg
// always passes idip and trigger; iso is never tested on the tag.

template <int NEta, int NCoeff, typename HIST_EFFMC>
class muon_insitu_efficiency_helper {
public:
  static_assert(NCoeff <= 4, "Chebyshev order > 3 not implemented");

  // flat index layout: [ idip | trigger | iso ]
  static constexpr int nID = NEta * 2 * NCoeff; // idip is charge split
  static constexpr int nHLT = NEta * NCoeff;    // trigger charge inclusive
  static constexpr int nIso = NEta * NCoeff;    // iso charge inclusive
  static constexpr int NSF = nID + nHLT + nIso;

  using tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<NSF>>;

  muon_insitu_efficiency_helper(HIST_EFFMC &&effMC, double ptmin, double ptmax,
                                double delta = 0.01, double effMC_max = 0.999)
      : effMC_(std::make_shared<const HIST_EFFMC>(std::move(effMC))),
        ptmin_(ptmin), ptmax_(ptmax), delta_(delta), effMC_max_(effMC_max) {}

  tensor_t operator()(float probe_pt, float probe_eta, int probe_charge,
                      float tag_pt, float tag_eta, int tag_charge, int category,
                      double nominal_weight) const {

    tensor_t res;
    res.setConstant(nominal_weight);

    std::array<double, NSF> grad;
    grad.fill(0.0);

    // tag leg always passes idip and trigger
    accumulate(grad, tag_pt, tag_eta, tag_charge, IDIP, true);
    accumulate(grad, tag_pt, tag_eta, tag_charge, TRIG, true);

    switch (category) {
    case 0: // nominal: probe passes idip, trigger, iso
      accumulate(grad, probe_pt, probe_eta, probe_charge, IDIP, true);
      accumulate(grad, probe_pt, probe_eta, probe_charge, TRIG, true);
      accumulate(grad, probe_pt, probe_eta, probe_charge, ISO, true);
      break;
    case 1: // failIso: probe passes idip, trigger, fails iso
      accumulate(grad, probe_pt, probe_eta, probe_charge, IDIP, true);
      accumulate(grad, probe_pt, probe_eta, probe_charge, TRIG, true);
      accumulate(grad, probe_pt, probe_eta, probe_charge, ISO, false);
      break;
    case 2: // failHLT: probe passes idip, fails trigger
      accumulate(grad, probe_pt, probe_eta, probe_charge, IDIP, true);
      accumulate(grad, probe_pt, probe_eta, probe_charge, TRIG, false);
      break;
    case 3: // failID: probe fails idip
      accumulate(grad, probe_pt, probe_eta, probe_charge, IDIP, false);
      break;
    default:
      break;
    }

    for (int i = 0; i < NSF; ++i) {
      if (grad[i] != 0.0) {
        double dlogk = delta_ * grad[i];
        // last-resort guard against exp overflow for pathological bins
        if (dlogk > 30.0)
          dlogk = 30.0;
        if (dlogk < -30.0)
          dlogk = -30.0;
        res(i) = nominal_weight * std::exp(dlogk);
      }
    }
    return res;
  }

private:
  enum Step { IDIP = 0, TRIG = 1, ISO = 2 };

  void accumulate(std::array<double, NSF> &grad, float pt, float eta,
                  int charge, Step step, bool pass) const {

    // Chebyshev argument on the probe pt window
    double p = pt;
    if (p < ptmin_)
      p = ptmin_;
    if (p > ptmax_)
      p = ptmax_;
    const double x = 2.0 * (p - ptmin_) / (ptmax_ - ptmin_) - 1.0;
    double T[4];
    T[0] = 1.0;
    T[1] = x;
    T[2] = 2.0 * x * x - 1.0;
    T[3] = 4.0 * x * x * x - 3.0 * x;

    // eta axis index (may be flow); the parameter decorrelation bin is the
    // same axis clamped into [0, NEta)
    const int eta_idx = effMC_->template axis<0>().index(eta);
    int b = eta_idx;
    if (b < 0)
      b = 0;
    if (b >= NEta)
      b = NEta - 1;

    const int qbit = (charge > 0) ? 1 : 0;

    double dlnf_ds;
    if (pass) {
      dlnf_ds = 1.0;
    } else {
      const int pt_idx = effMC_->template axis<1>().index(pt);
      const int charge_idx = effMC_->template axis<2>().index(charge);
      double e =
          effMC_->at(eta_idx, pt_idx, charge_idx, stepIndex(step)).value();
      if (!(e > 0.0)) {
        return; // empty effMC -> no variation for this step
      }
      if (e > effMC_max_)
        e = effMC_max_; // bound the fail-leg derivative for eMC -> 1
      dlnf_ds = -e / (1.0 - e);
    }

    int base;
    int blockBin;
    if (step == IDIP) {
      base = 0;
      blockBin = qbit * NEta + b; // idip is charge split
    } else if (step == TRIG) {
      base = nID;
      blockBin = b;
    } else {
      base = nID + nHLT;
      blockBin = b;
    }
    const int off = base + blockBin * NCoeff;
    for (int k = 0; k < NCoeff; ++k) {
      grad[off + k] += T[k] * dlnf_ds;
    }
  }

  int stepIndex(Step step) const {
    if (step == IDIP)
      return idx_idip_;
    if (step == TRIG)
      return idx_trig_;
    return idx_iso_;
  }

  std::shared_ptr<const HIST_EFFMC> effMC_;
  double ptmin_;
  double ptmax_;
  double delta_;
  double effMC_max_;
  // cache the slow string-category lookups
  int idx_idip_ = effMC_->template axis<3>().index("idip");
  int idx_trig_ = effMC_->template axis<3>().index("trigger");
  int idx_iso_ = effMC_->template axis<3>().index("iso");
};

} // namespace wrem

#endif
