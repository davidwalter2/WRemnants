#ifndef WREMNANTS_MUON_EFFICIENCIES_INSITU_H
#define WREMNANTS_MUON_EFFICIENCIES_INSITU_H

#include <array>
#include <boost/histogram/axis.hpp>
#include <cmath>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <memory>

namespace wrem {

// In-situ muon efficiency helper for the W/Z muon analyses.
//
// Floats the ID (idip), HLT (trigger) and Iso efficiency data/MC scale
// factors as Chebyshev polynomials, decorrelated per NEta-bin probe-eta:
//   IDIP    : 1D in pt, charge-dependent  -> NCoeffPt per (eta, q)
//   Trigger : 2D in (pt, ut), charge-dep  -> NCoeffPt*NCoeffUt per (eta, q)
//   Iso     : 2D in (pt, ut), charge-incl -> NCoeffPt*NCoeffUt per (eta)
// matching the SMP-23-002 measurement variables.
//
//   xtilde_pt = 2*(clamp(pt,ptmin,ptmax)-ptmin)/(ptmax-ptmin) - 1
//   xtilde_ut = 2*(clamp(ut,utmin,utmax)-utmin)/(utmax-utmin) - 1
//   SF_X = exp( sum_k [, m] theta_{X,etaBin,[q],k[,m]} T_k(x_pt) [T_m(x_ut)] )
// A leg PASSING step X contributes a factor SF_X; FAILING contributes
//   (1 - SF_X * eMC) / (1 - eMC),   eMC = effMC_X(pt, eta, q[, ut]).
// At theta = 0 every factor = 1 -> nominal unchanged.
//
// rabbit treats every output-tensor bin as an independent linearised
// (log-normal) nuisance around the nominal. We store per coefficient the
// analytic first-order response:
//   res(c) = w_nom * exp( delta * d lnW / d theta_c )
// with d ln f / d s = +1 (pass) or -eMC/(1-eMC) (fail). Coefficients are
// unconstrained, so delta only sets nuisance units (theta_c = delta * n_c).
// delta is kept small because the fail-leg derivative is large for high
// efficiency (idip eMC ~ 0.99). effMC clamped to [0, effMC_max].
//
// category: 0 = nominal, 1 = failIso, 2 = failHLT, 3 = failID.
//
// Two concrete helpers share the same implementation through a common base:
//   - muon_insitu_efficiency_helper          (two-leg, Z dilepton): probe +
//       tag muon variables; the tag leg always passes idip & trigger (the
//       dilepton tag-and-probe topology), iso is never tested on the tag.
//   - muon_insitu_efficiency_helper_singleleg (single-leg, W single muon):
//       only the probe (the single good muon) contributes; there is no tag
//       leg. Use category 0 (nominal) for signal muons that pass idip,
//       trigger and iso.
// Each derived class exposes a single operator() (distinct signatures) so the
// RDF backend dispatches unambiguously.

template <int NEta, int NCoeffPt, int NCoeffUt, typename HIST_IDIP,
          typename HIST_TRIG, typename HIST_ISO>
class muon_insitu_efficiency_helper_base {
public:
  static_assert(NCoeffPt <= 4, "Chebyshev pt order > 3 not implemented");
  static_assert(NCoeffUt <= 4, "Chebyshev ut order > 3 not implemented");

  // flat index layout: [ idip | trigger | iso ]
  static constexpr int NCoeff2D = NCoeffPt * NCoeffUt;
  static constexpr int nID = NEta * 2 * NCoeffPt;  // (eta, q) x pt
  static constexpr int nHLT = NEta * 2 * NCoeff2D; // (eta, q) x pt x ut
  static constexpr int nIso = NEta * NCoeff2D;     // eta x pt x ut
  static constexpr int NSF = nID + nHLT + nIso;

  using tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<NSF>>;

  muon_insitu_efficiency_helper_base(HIST_IDIP &&effMC_idip,
                                     HIST_TRIG &&effMC_trig,
                                     HIST_ISO &&effMC_iso, double ptmin,
                                     double ptmax, double utmin, double utmax,
                                     double delta = 0.01,
                                     double effMC_max = 0.999)
      : effMC_idip_(std::make_shared<const HIST_IDIP>(std::move(effMC_idip))),
        effMC_trig_(std::make_shared<const HIST_TRIG>(std::move(effMC_trig))),
        effMC_iso_(std::make_shared<const HIST_ISO>(std::move(effMC_iso))),
        ptmin_(ptmin), ptmax_(ptmax), utmin_(utmin), utmax_(utmax),
        delta_(delta), effMC_max_(effMC_max) {}

protected:
  enum Step { IDIP = 0, TRIG = 1, ISO = 2 };

  // Shared response computation. WithTag adds the (always-passing idip &
  // trigger) tag leg at compile time; for the single-leg helper it is elided.
  template <bool WithTag>
  tensor_t compute(float probe_pt, float probe_eta, int probe_charge,
                   float probe_ut, float tag_pt, float tag_eta, int tag_charge,
                   float tag_ut, int category, double nominal_weight) const {

    tensor_t res;
    res.setConstant(nominal_weight);

    std::array<double, NSF> grad;
    grad.fill(0.0);

    // tag leg always passes idip and trigger (two-leg topology only)
    if constexpr (WithTag) {
      accumulate(grad, tag_pt, tag_eta, tag_charge, tag_ut, IDIP, true);
      accumulate(grad, tag_pt, tag_eta, tag_charge, tag_ut, TRIG, true);
    }

    switch (category) {
    case 0: // nominal: probe passes idip, trigger, iso
      accumulate(grad, probe_pt, probe_eta, probe_charge, probe_ut, IDIP, true);
      accumulate(grad, probe_pt, probe_eta, probe_charge, probe_ut, TRIG, true);
      accumulate(grad, probe_pt, probe_eta, probe_charge, probe_ut, ISO, true);
      break;
    case 1: // failIso: probe passes idip, trigger, fails iso
      accumulate(grad, probe_pt, probe_eta, probe_charge, probe_ut, IDIP, true);
      accumulate(grad, probe_pt, probe_eta, probe_charge, probe_ut, TRIG, true);
      accumulate(grad, probe_pt, probe_eta, probe_charge, probe_ut, ISO, false);
      break;
    case 2: // failHLT: probe passes idip, fails trigger
      accumulate(grad, probe_pt, probe_eta, probe_charge, probe_ut, IDIP, true);
      accumulate(grad, probe_pt, probe_eta, probe_charge, probe_ut, TRIG,
                 false);
      break;
    case 3: // failID: probe fails idip
      accumulate(grad, probe_pt, probe_eta, probe_charge, probe_ut, IDIP,
                 false);
      break;
    default:
      break;
    }

    for (int i = 0; i < NSF; ++i) {
      if (grad[i] != 0.0) {
        double dlogk = delta_ * grad[i];
        if (dlogk > 30.0)
          dlogk = 30.0;
        if (dlogk < -30.0)
          dlogk = -30.0;
        res(i) = nominal_weight * std::exp(dlogk);
      }
    }
    return res;
  }

  void accumulate(std::array<double, NSF> &grad, float pt, float eta,
                  int charge, float ut, Step step, bool pass) const {

    // Chebyshev arguments
    double p = pt;
    if (p < ptmin_)
      p = ptmin_;
    if (p > ptmax_)
      p = ptmax_;
    const double xpt = 2.0 * (p - ptmin_) / (ptmax_ - ptmin_) - 1.0;
    double Tpt[4];
    Tpt[0] = 1.0;
    Tpt[1] = xpt;
    Tpt[2] = 2.0 * xpt * xpt - 1.0;
    Tpt[3] = 4.0 * xpt * xpt * xpt - 3.0 * xpt;

    double Tut[4] = {1.0, 0.0, 0.0, 0.0};
    if (step != IDIP) {
      double u = ut;
      if (u < utmin_)
        u = utmin_;
      if (u > utmax_)
        u = utmax_;
      const double xut = 2.0 * (u - utmin_) / (utmax_ - utmin_) - 1.0;
      Tut[0] = 1.0;
      Tut[1] = xut;
      Tut[2] = 2.0 * xut * xut - 1.0;
      Tut[3] = 4.0 * xut * xut * xut - 3.0 * xut;
    }

    // eta decorrelation bin (clamped, NEta-grid)
    const int eta_idx_idip = effMC_idip_->template axis<0>().index(eta);
    int b = eta_idx_idip;
    if (b < 0)
      b = 0;
    if (b >= NEta)
      b = NEta - 1;

    const int qbit = (charge > 0) ? 1 : 0;

    // effMC lookup per step
    double dlnf_ds;
    if (pass) {
      dlnf_ds = 1.0;
    } else {
      double e = 0.0;
      if (step == IDIP) {
        const int eta_i = effMC_idip_->template axis<0>().index(eta);
        const int pt_i = effMC_idip_->template axis<1>().index(pt);
        const int q_i = effMC_idip_->template axis<2>().index(charge);
        e = effMC_idip_->at(eta_i, pt_i, q_i).value();
      } else if (step == TRIG) {
        const int eta_i = effMC_trig_->template axis<0>().index(eta);
        const int pt_i = effMC_trig_->template axis<1>().index(pt);
        const int q_i = effMC_trig_->template axis<2>().index(charge);
        const int ut_i = effMC_trig_->template axis<3>().index(ut);
        e = effMC_trig_->at(eta_i, pt_i, q_i, ut_i).value();
      } else { // ISO (charge-inclusive)
        const int eta_i = effMC_iso_->template axis<0>().index(eta);
        const int pt_i = effMC_iso_->template axis<1>().index(pt);
        const int ut_i = effMC_iso_->template axis<2>().index(ut);
        e = effMC_iso_->at(eta_i, pt_i, ut_i).value();
      }
      if (!(e > 0.0)) {
        return; // empty effMC -> no variation for this step
      }
      if (e > effMC_max_)
        e = effMC_max_;
      dlnf_ds = -e / (1.0 - e);
    }

    // index layout per step block
    if (step == IDIP) {
      // (q, eta) x pt
      const int blockBin = qbit * NEta + b;
      const int off = 0 + blockBin * NCoeffPt;
      for (int k = 0; k < NCoeffPt; ++k) {
        grad[off + k] += Tpt[k] * dlnf_ds;
      }
    } else if (step == TRIG) {
      // (q, eta) x pt x ut
      const int blockBin = qbit * NEta + b;
      const int off = nID + blockBin * NCoeff2D;
      for (int k = 0; k < NCoeffPt; ++k) {
        for (int m = 0; m < NCoeffUt; ++m) {
          grad[off + k * NCoeffUt + m] += Tpt[k] * Tut[m] * dlnf_ds;
        }
      }
    } else { // ISO
      // eta x pt x ut (charge-inclusive)
      const int blockBin = b;
      const int off = nID + nHLT + blockBin * NCoeff2D;
      for (int k = 0; k < NCoeffPt; ++k) {
        for (int m = 0; m < NCoeffUt; ++m) {
          grad[off + k * NCoeffUt + m] += Tpt[k] * Tut[m] * dlnf_ds;
        }
      }
    }
  }

  std::shared_ptr<const HIST_IDIP> effMC_idip_;
  std::shared_ptr<const HIST_TRIG> effMC_trig_;
  std::shared_ptr<const HIST_ISO> effMC_iso_;
  double ptmin_;
  double ptmax_;
  double utmin_;
  double utmax_;
  double delta_;
  double effMC_max_;
};

// Two-leg helper (Z dilepton tag-and-probe).
template <int NEta, int NCoeffPt, int NCoeffUt, typename HIST_IDIP,
          typename HIST_TRIG, typename HIST_ISO>
class muon_insitu_efficiency_helper
    : public muon_insitu_efficiency_helper_base<
          NEta, NCoeffPt, NCoeffUt, HIST_IDIP, HIST_TRIG, HIST_ISO> {
public:
  using Base =
      muon_insitu_efficiency_helper_base<NEta, NCoeffPt, NCoeffUt, HIST_IDIP,
                                         HIST_TRIG, HIST_ISO>;
  using Base::Base;
  using typename Base::tensor_t;

  tensor_t operator()(float probe_pt, float probe_eta, int probe_charge,
                      float probe_ut, float tag_pt, float tag_eta,
                      int tag_charge, float tag_ut, int category,
                      double nominal_weight) const {
    return this->template compute<true>(probe_pt, probe_eta, probe_charge,
                                        probe_ut, tag_pt, tag_eta, tag_charge,
                                        tag_ut, category, nominal_weight);
  }
};

// Single-leg helper (W single muon): no tag leg.
template <int NEta, int NCoeffPt, int NCoeffUt, typename HIST_IDIP,
          typename HIST_TRIG, typename HIST_ISO>
class muon_insitu_efficiency_helper_singleleg
    : public muon_insitu_efficiency_helper_base<
          NEta, NCoeffPt, NCoeffUt, HIST_IDIP, HIST_TRIG, HIST_ISO> {
public:
  using Base =
      muon_insitu_efficiency_helper_base<NEta, NCoeffPt, NCoeffUt, HIST_IDIP,
                                         HIST_TRIG, HIST_ISO>;
  using Base::Base;
  using typename Base::tensor_t;

  tensor_t operator()(float probe_pt, float probe_eta, int probe_charge,
                      float probe_ut, int category,
                      double nominal_weight) const {
    return this->template compute<false>(probe_pt, probe_eta, probe_charge,
                                         probe_ut, 0.0f, 0.0f, 0, 0.0f,
                                         category, nominal_weight);
  }
};

} // namespace wrem

#endif
