#ifndef WREMNANTS_MUON_EFFICIENCIES_INSITU_H
#define WREMNANTS_MUON_EFFICIENCIES_INSITU_H

#include <array>
#include <boost/histogram/axis.hpp>
#include <cmath>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <memory>
#include <vector>

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
//   s_X = sum_k [, m] theta_{X,etaBin,[q],k[,m]} T_k(x_pt) [T_m(x_ut)]
//
// LOGIT parameterisation of the data efficiency (bounded in (0,1), so the
// fail factor never hits the (1 - SF*eMC) -> 0 singularity): with the MC
// efficiency eMC = effMC_X(pt,eta,q[,ut]) (clamped to (0, effMC_max] so that
// 1-eMC > 0), expS = exp(s_X), and the common denominator
//   D = eMC*expS + (1 - eMC),
// the data efficiency is  eps = eMC*expS / D = sigmoid(s + logit(eMC)) in
// (0,1). A leg PASSING step X contributes  f_pass = expS/D ( = eps/eMC );
// FAILING contributes  f_fail = 1/D ( = (1-eps)/(1-eMC) ). At theta = 0:
// expS=1, D=1, every factor = 1 -> nominal unchanged.
//
// rabbit treats every output-tensor bin as an independent linearised
// (log-normal) nuisance around the nominal. We store per coefficient the
// analytic first-order response:
//   res(c) = w_nom * exp( delta * d lnW / d theta_c )
// with d ln f / d s = (1 - eps) (pass) or -eps (fail) -- both in [-1,1], so the
// linearisation is well-conditioned and NO clamp on the gradient is needed (the
// old exp(s) form had a divergent fail derivative -SF*eMC/(1-SF*eMC)).
// Coefficients are unconstrained, so delta only sets nuisance units
// (theta_c = delta * n_c). effMC is clamped to (0, effMC_max].
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

  // theta_central (optional, last arg): the accumulated central Chebyshev
  // coefficients for the iterative fit. Empty -> all zero (SF=1, linearise
  // around the nominal == iteration 0). Otherwise the variations are
  // linearised, and the MC reweighted (via central_weight), around
  // SF = exp(s(theta_central)).
  muon_insitu_efficiency_helper_base(
      HIST_IDIP &&effMC_idip, HIST_TRIG &&effMC_trig, HIST_ISO &&effMC_iso,
      double ptmin, double ptmax, double utmin, double utmax,
      double delta = 0.01, double effMC_max = 0.9999,
      const std::vector<double> &theta_central = {})
      : effMC_idip_(std::make_shared<const HIST_IDIP>(std::move(effMC_idip))),
        effMC_trig_(std::make_shared<const HIST_TRIG>(std::move(effMC_trig))),
        effMC_iso_(std::make_shared<const HIST_ISO>(std::move(effMC_iso))),
        ptmin_(ptmin), ptmax_(ptmax), utmin_(utmin), utmax_(utmax),
        delta_(delta), effMC_max_(effMC_max) {
    if (theta_central.empty()) {
      theta_central_.fill(0.0);
    } else {
      for (int i = 0; i < NSF; ++i)
        theta_central_[i] = theta_central[i];
    }
  }

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

    // delta*grad is bounded by construction (|d ln f/ds| <= 1, |T| <= 1), so no
    // overflow guard on the exponent is needed.
    for (int i = 0; i < NSF; ++i) {
      if (grad[i] != 0.0)
        res(i) = nominal_weight * std::exp(delta_ * grad[i]);
    }
    return res;
  }

  // Compute the Chebyshev factors Tpt/Tut, the flat block offset, and
  // expS = exp(clamp(s(theta_central),±30)) for one leg/step (s is the log of
  // the SF in the old basis; here it is the logit shift). Shared by accumulate
  // (gradient) and central_factor (reweight) so both linearise around the SAME
  // point -> self-consistent fit.
  void eval_leg(float pt, float eta, int charge, float ut, Step step,
                double Tpt[4], double Tut[4], int &off, double &expS) const {
    double p = pt;
    if (p < ptmin_)
      p = ptmin_;
    if (p > ptmax_)
      p = ptmax_;
    const double xpt = 2.0 * (p - ptmin_) / (ptmax_ - ptmin_) - 1.0;
    Tpt[0] = 1.0;
    Tpt[1] = xpt;
    Tpt[2] = 2.0 * xpt * xpt - 1.0;
    Tpt[3] = 4.0 * xpt * xpt * xpt - 3.0 * xpt;

    Tut[0] = 1.0;
    Tut[1] = 0.0;
    Tut[2] = 0.0;
    Tut[3] = 0.0;
    if (step != IDIP) {
      double u = ut;
      if (u < utmin_)
        u = utmin_;
      if (u > utmax_)
        u = utmax_;
      const double xut = 2.0 * (u - utmin_) / (utmax_ - utmin_) - 1.0;
      Tut[1] = xut;
      Tut[2] = 2.0 * xut * xut - 1.0;
      Tut[3] = 4.0 * xut * xut * xut - 3.0 * xut;
    }

    int b = effMC_idip_->template axis<0>().index(eta);
    if (b < 0)
      b = 0;
    if (b >= NEta)
      b = NEta - 1;
    const int qbit = (charge > 0) ? 1 : 0;

    double s = 0.0;
    if (step == IDIP) {
      off = (qbit * NEta + b) * NCoeffPt;
      for (int k = 0; k < NCoeffPt; ++k)
        s += theta_central_[off + k] * Tpt[k];
    } else if (step == TRIG) {
      off = nID + (qbit * NEta + b) * NCoeff2D;
      for (int k = 0; k < NCoeffPt; ++k)
        for (int m = 0; m < NCoeffUt; ++m)
          s += theta_central_[off + k * NCoeffUt + m] * Tpt[k] * Tut[m];
    } else { // ISO (charge-inclusive)
      off = nID + nHLT + b * NCoeff2D;
      for (int k = 0; k < NCoeffPt; ++k)
        for (int m = 0; m < NCoeffUt; ++m)
          s += theta_central_[off + k * NCoeffUt + m] * Tpt[k] * Tut[m];
    }
    if (s > 30.0)
      s = 30.0;
    if (s < -30.0)
      s = -30.0;
    expS = std::exp(s);
  }

  // effMC lookup for a fail leg; returns e clamped to (0, effMC_max_], or 0.0
  // to signal "empty effMC -> skip this step".
  double lookup_effMC(Step step, float pt, float eta, int charge,
                      float ut) const {
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
    if (!(e > 0.0))
      return 0.0;
    if (e > effMC_max_)
      e = effMC_max_;
    return e;
  }

  // Accumulate the linearised gradient dlnW/dtheta_c at theta_central for one
  // leg/step. Under the logit parameterisation both derivatives are bounded:
  //   d ln f / d s = (1 - eps) (pass) or -eps (fail),
  // with eps = eMC*expS/(eMC*expS + 1-eMC) the data efficiency in (0,1). At
  // theta_central=0 (expS=1) -> (1-eMC) (pass) / -eMC (fail). Both legs now
  // need effMC (the old pass derivative +1 was effMC-independent).
  void accumulate(std::array<double, NSF> &grad, float pt, float eta,
                  int charge, float ut, Step step, bool pass) const {
    double Tpt[4], Tut[4], expS;
    int off;
    eval_leg(pt, eta, charge, ut, step, Tpt, Tut, off, expS);

    double e = lookup_effMC(step, pt, eta, charge, ut);
    if (e <= 0.0)
      return; // empty effMC -> no variation for this step
    const double eps = e * expS / (e * expS + (1.0 - e));
    const double dlnf_ds = pass ? (1.0 - eps) : -eps;

    if (step == IDIP) {
      for (int k = 0; k < NCoeffPt; ++k)
        grad[off + k] += Tpt[k] * dlnf_ds;
    } else {
      for (int k = 0; k < NCoeffPt; ++k)
        for (int m = 0; m < NCoeffUt; ++m)
          grad[off + k * NCoeffUt + m] += Tpt[k] * Tut[m] * dlnf_ds;
    }
  }

  // Central per-leg/step factor f_X(theta_central): logit closed forms
  // f_pass = expS/D, f_fail = 1/D with D = e*expS + (1-e). Equals 1 at
  // theta_central = 0. The closed forms are used directly (not (1-eps)/(1-e))
  // so e -> 1 is finite (the fail category is empty in MC there anyway).
  double central_factor(float pt, float eta, int charge, float ut, Step step,
                        bool pass) const {
    double Tpt[4], Tut[4], expS;
    int off;
    eval_leg(pt, eta, charge, ut, step, Tpt, Tut, off, expS);
    double e = lookup_effMC(step, pt, eta, charge, ut);
    if (e <= 0.0)
      return 1.0; // empty effMC -> no reweight
    const double D = e * expS + (1.0 - e);
    return pass ? (expS / D) : (1.0 / D);
  }

  // Central MC reweight W(theta_central) = product of per-leg/step factors over
  // the tested steps of the category. 1.0 at theta_central = 0.
  template <bool WithTag>
  double central_weight(float probe_pt, float probe_eta, int probe_charge,
                        float probe_ut, float tag_pt, float tag_eta,
                        int tag_charge, float tag_ut, int category) const {
    double W = 1.0;
    if constexpr (WithTag) {
      W *= central_factor(tag_pt, tag_eta, tag_charge, tag_ut, IDIP, true);
      W *= central_factor(tag_pt, tag_eta, tag_charge, tag_ut, TRIG, true);
    }
    switch (category) {
    case 0: // nominal
      W *= central_factor(probe_pt, probe_eta, probe_charge, probe_ut, IDIP,
                          true);
      W *= central_factor(probe_pt, probe_eta, probe_charge, probe_ut, TRIG,
                          true);
      W *= central_factor(probe_pt, probe_eta, probe_charge, probe_ut, ISO,
                          true);
      break;
    case 1: // failIso
      W *= central_factor(probe_pt, probe_eta, probe_charge, probe_ut, IDIP,
                          true);
      W *= central_factor(probe_pt, probe_eta, probe_charge, probe_ut, TRIG,
                          true);
      W *= central_factor(probe_pt, probe_eta, probe_charge, probe_ut, ISO,
                          false);
      break;
    case 2: // failHLT
      W *= central_factor(probe_pt, probe_eta, probe_charge, probe_ut, IDIP,
                          true);
      W *= central_factor(probe_pt, probe_eta, probe_charge, probe_ut, TRIG,
                          false);
      break;
    case 3: // failID
      W *= central_factor(probe_pt, probe_eta, probe_charge, probe_ut, IDIP,
                          false);
      break;
    default:
      break;
    }
    // W is a product of bounded positive factors (f_pass in (0,1/e],
    // f_fail in (0,1/(1-e)], with expS bounded by the +-30 clamp in eval_leg),
    // so it stays finite and positive without an explicit bound.
    return W;
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
  std::array<double, NSF> theta_central_;
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

// Central reweight functor (two-leg): returns the scalar W(theta_central) to
// fold into nominal_weight before the templates/tensors are filled. Shares the
// base (effMC + theta_central) with the tensor helper.
template <int NEta, int NCoeffPt, int NCoeffUt, typename HIST_IDIP,
          typename HIST_TRIG, typename HIST_ISO>
class muon_insitu_central_weight_helper
    : public muon_insitu_efficiency_helper_base<
          NEta, NCoeffPt, NCoeffUt, HIST_IDIP, HIST_TRIG, HIST_ISO> {
public:
  using Base =
      muon_insitu_efficiency_helper_base<NEta, NCoeffPt, NCoeffUt, HIST_IDIP,
                                         HIST_TRIG, HIST_ISO>;
  using Base::Base;

  double operator()(float probe_pt, float probe_eta, int probe_charge,
                    float probe_ut, float tag_pt, float tag_eta, int tag_charge,
                    float tag_ut, int category) const {
    return this->template central_weight<true>(
        probe_pt, probe_eta, probe_charge, probe_ut, tag_pt, tag_eta,
        tag_charge, tag_ut, category);
  }
};

// Central reweight functor (single-leg): no tag leg.
template <int NEta, int NCoeffPt, int NCoeffUt, typename HIST_IDIP,
          typename HIST_TRIG, typename HIST_ISO>
class muon_insitu_central_weight_helper_singleleg
    : public muon_insitu_efficiency_helper_base<
          NEta, NCoeffPt, NCoeffUt, HIST_IDIP, HIST_TRIG, HIST_ISO> {
public:
  using Base =
      muon_insitu_efficiency_helper_base<NEta, NCoeffPt, NCoeffUt, HIST_IDIP,
                                         HIST_TRIG, HIST_ISO>;
  using Base::Base;

  double operator()(float probe_pt, float probe_eta, int probe_charge,
                    float probe_ut, int category) const {
    return this->template central_weight<false>(probe_pt, probe_eta,
                                                probe_charge, probe_ut, 0.0f,
                                                0.0f, 0, 0.0f, category);
  }
};

} // namespace wrem

#endif
