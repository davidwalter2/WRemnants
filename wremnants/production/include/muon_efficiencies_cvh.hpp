#ifndef WREMNANTS_MUON_EFFICIENCIES_CVH_H
#define WREMNANTS_MUON_EFFICIENCIES_CVH_H

#include <ROOT/RVec.hxx>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>

// Treatment of the CVH-refit efficiency holes.
//
// Two badly aligned modules in the real data tracker alignment kill the tracks
// that cross them, each in its own way:
//
//   * the glued (double-sided) module in TIB layer 2, detId 369141860, whose
//     mono face is rotated ~0.75 rad from ideal, corrupting the composite frame
//     that anchors the CVH fit (fixed in WMass/cmssw PR#46, commit eb96caef);
//
//   * the r-phi sensor on TID disk -2, detId 402666798, sitting 5.3 mm out of
//     its disk plane, which makes the CVH propagation abort mid-track (still
//     uncorrected).
//
// The effect is DATA-ONLY -- MC is reconstructed with ideal geometry and has no
// such pathology.
//
// Both affected regions are localized rectangles in (eta, phi'), measured
// directly from the single-muon CVH refit efficiency in Z->mumu events with
//   mz_dilepton.py --cvhEfficiencyHists,
// and written out by
// scripts/analysisTools/w_mass_13TeV/make_cvh_efficiency_sf.py. Everywhere
// outside them the data/MC agreement is at the 1e-5 level.
//
// phi' is the azimuth at the module rather than at the vertex: a muon of
// charge q and transverse momentum pt is bent by q*C/pt on its way out to the
// module radius r, with C = 0.3*B*r/2. Undoing it (see cvh_module_phi) makes
// the regions charge- and pt-independent, and is what keeps the sharp edges of
// the module from being smeared over the two charges.
//
// Two treatments are provided:
//
//  * cvh_muon_sf (--cvhBadModules sf, the histmaker default): the measured
//    SF = eps_data/eps_MC, applied as a per-muon weight on MC only. It has to
//    be paired with cvh_antiveto_helper below, which carries the mirror effect
//    on the muons that are meant *not* to be reconstructed -- currently done in
//    mw_with_mu_eta_pt.py only, the dilepton selections having no anti-veto SF
//    to correct.
//
//  * cvh_muon_in_bad_module (--cvhBadModules veto): a geometric veto. Events
//    with any muon crossing one of the regions are dropped from DATA AND MC
//    alike, so the affected phase space simply leaves the measurement and no
//    efficiency correction is needed at all. Because the decision is made from
//    the uncorrected (standard tracker fit) kinematics it is blind to whether
//    the CVH refit succeeded, hence identical in data and MC by construction.
//    Costs ~1.3% of the muons, and is the cross-check on the SF treatment.

namespace wrem {

// Bending to the module radius, C = 0.3*B*r/2 in rad*GeV, i.e. C = 0.57*r with
// r in metres for B = 3.8 T. It belongs to the module, not to the analysis, so
// each region below carries its own value (bendC), taken from the radius at
// which its module sits: rho = 34.96 cm for the TIB one (C = 0.199) and
// 36.97 cm for the TID one (C = 0.211). make_cvh_efficiency_sf.py --fitBending
// measures the same quantity from the charge splitting of each hole and agrees
// with both within ~1 sigma, but only usefully constrains the deep TIB one.
//
// This constant is just the fallback for callers that ask for a module phi
// without naming a region; the regions themselves never use it.
inline constexpr double cvh_bending_C = 0.199;

// Azimuth of the muon where it crosses the module, in [0, 2*pi).
// Accepts phi either in [-pi,pi] (nanoAOD) or in [0,2*pi).
inline float cvh_module_phi(float phi, int charge, float pt,
                            double bendC = cvh_bending_C) {
  constexpr float twopi = 2.f * static_cast<float>(M_PI);
  float p = phi - static_cast<float>(charge * bendC) / pt;
  while (p < 0.f)
    p += twopi;
  while (p >= twopi)
    p -= twopi;
  return p;
}

// A rectangular region of the (eta, phi') plane holding the measured SF map,
// row-major in [ieta][iphi]. phi' is built with this region's own bendC.
struct CvhSFRegion {
  double etaLow, etaBinWidth;
  int nEta;
  double phiLow, phiBinWidth;
  int nPhi;
  double bendC;
  const double *sf;
};

// The same regions as plain (eta, phi') rectangles, for the geometric veto.
struct CvhModuleRegion {
  double etaLow, etaHigh, phiLow, phiHigh, bendC;
};

// clang-format off
// BEGIN GENERATED -- do not edit by hand, see make_cvh_efficiency_sf.py

// hotspot1 (TIB-L2, detId 369141860)
//   eta in [0.00, 0.75), phi' in [0.698, 1.047), C = 0.199 (r = 35 cm)
inline const double cvh_sf_region0[15 * 4] = {
    1.000, 0.991, 0.987, 1.000, // eta ~ +0.025
    1.000, 0.973, 0.975, 1.000, // eta ~ +0.075
    1.000, 0.933, 0.940, 1.000, // eta ~ +0.125
    1.000, 0.892, 0.886, 0.995, // eta ~ +0.175
    0.995, 0.800, 0.770, 0.990, // eta ~ +0.225
    0.993, 0.744, 0.635, 0.978, // eta ~ +0.275
    0.993, 0.690, 0.503, 0.960, // eta ~ +0.325
    0.993, 0.717, 0.419, 0.943, // eta ~ +0.375
    0.995, 0.771, 0.394, 0.928, // eta ~ +0.425
    0.997, 0.860, 0.490, 0.926, // eta ~ +0.475
    1.000, 0.925, 0.635, 0.931, // eta ~ +0.525
    1.000, 0.971, 0.782, 0.950, // eta ~ +0.575
    1.000, 0.991, 0.912, 0.967, // eta ~ +0.625
    1.000, 1.000, 0.966, 0.989, // eta ~ +0.675
    1.000, 1.000, 0.993, 0.994, // eta ~ +0.725
};

// hotspot2 (TID-2, detId 402666798)
//   eta in [-1.80, -1.45), phi' in [5.061, 5.411), C = 0.211 (r = 37 cm)
inline const double cvh_sf_region1[7 * 4] = {
    1.000, 1.000, 1.000, 0.983, // eta ~ -1.775
    0.992, 0.990, 0.990, 0.922, // eta ~ -1.725
    0.986, 0.968, 0.964, 0.872, // eta ~ -1.675
    0.958, 0.953, 0.939, 0.887, // eta ~ -1.625
    0.931, 0.923, 0.885, 0.902, // eta ~ -1.575
    0.965, 0.962, 0.943, 0.956, // eta ~ -1.525
    1.000, 1.000, 0.991, 0.994, // eta ~ -1.475
};

// last column of each row is bendC = 0.3*B*r/2, the bending to that module
inline const CvhSFRegion cvh_sf_regions[] = {
    {0.0000, 0.0500, 15, 0.6981, 0.0873, 4, 0.199, cvh_sf_region0},
    {-1.8000, 0.0500, 7, 5.0615, 0.0873, 4, 0.211, cvh_sf_region1},
};

// Bounding rectangles of the same regions, used by the geometric veto.
inline const CvhModuleRegion cvh_bad_modules[] = {
    {0.0000, 0.7500, 0.6981, 1.0472, 0.199}, // hotspot1 (TIB-L2, detId 369141860)
    {-1.8000, -1.4500, 5.0615, 5.4105, 0.211}, // hotspot2 (TID-2, detId 402666798)
};

// END GENERATED
// clang-format on

// True if the muon crosses one of the pathological modules. Meant to be fed
// the UNCORRECTED (standard tracker fit) kinematics, which exist whether or
// not the CVH refit succeeded; that is what makes the decision identical in
// data and MC. phi may be given in [-pi,pi] or [0,2*pi).
inline bool cvh_muon_in_bad_module(float eta, float phi, int charge, float pt) {
  for (const auto &r : cvh_bad_modules) {
    if (eta < r.etaLow || eta >= r.etaHigh)
      continue; // cheap test first, phi' costs a division per region
    const float phiModule = cvh_module_phi(phi, charge, pt, r.bendC);
    if (phiModule >= r.phiLow && phiModule < r.phiHigh)
      return true;
  }
  return false;
}

// Per-muon mask over a muon collection, for use in an event-level veto.
// Templated on the charge container only, whose element type varies with the
// nanoAOD version (Muon_charge is not necessarily an RVec<int>).
template <typename Vcharge>
inline ROOT::VecOps::RVec<bool> cvh_muons_in_bad_module(
    const ROOT::VecOps::RVec<float> &eta, const ROOT::VecOps::RVec<float> &phi,
    const Vcharge &charge, const ROOT::VecOps::RVec<float> &pt) {
  ROOT::VecOps::RVec<bool> out(eta.size(), false);
  for (std::size_t i = 0; i < eta.size(); ++i)
    out[i] = cvh_muon_in_bad_module(eta[i], phi[i], static_cast<int>(charge[i]),
                                    pt[i]);
  return out;
}

// Per-muon SF: eps_data/eps_MC, = 1 outside the affected regions. Apply to MC
// only, as a multiplicative weight, once per reconstructed muon. Alternative
// to the veto above, not to be combined with it.
inline double cvh_muon_sf(float eta, float phi, int charge, float pt) {
  for (const auto &r : cvh_sf_regions) {
    // floor, not truncation: the regions can sit at negative eta
    const int ieta =
        static_cast<int>(std::floor((eta - r.etaLow) / r.etaBinWidth));
    if (ieta < 0 || ieta >= r.nEta)
      continue;
    // each module is unbent to its own radius, so phi' is per region
    const float phiModule = cvh_module_phi(phi, charge, pt, r.bendC);
    const int iphi =
        static_cast<int>(std::floor((phiModule - r.phiLow) / r.phiBinWidth));
    if (iphi < 0 || iphi >= r.nPhi)
      continue;
    return r.sf[ieta * r.nPhi + iphi];
  }
  return 1.0;
}

// Event-level weight for a dimuon (Z) final state: both reconstructed muons
// must survive, so the per-muon SFs multiply.
inline double cvh_dimuon_sf(float eta1, float phi1, int charge1, float pt1,
                            float eta2, float phi2, int charge2, float pt2) {
  return cvh_muon_sf(eta1, phi1, charge1, pt1) *
         cvh_muon_sf(eta2, phi2, charge2, pt2);
}

// ---------------------------------------------------------------------------
// The same holes seen from the other side: the anti-veto.
//
// A dimuon event enters the single-muon (W) selection when the second muon
// escapes the veto. In data that can also happen because its CVH refit failed
// -- a muon with charge = -99 is not in 'vetoMuons' -- while in MC, with ideal
// geometry, it cannot. So DY leaks into the W selection in data only.
//
// The veto SF machinery already reweights those events by the anti-veto SF
//   a = (1 - eps_MC^veto * SF) / (1 - eps_MC^veto) = (1 - eps_d) / (1 - eps_MC)
// evaluated on the postFSR gen muon that was not reconstructed. Folding the
// refit efficiency in means replacing eps_d by eps_d * s, with s = cvh_muon_sf
// the measured refit SF (eps_MC^cvh = 1 for ideal geometry):
//
//   a' = (1 - eps_d * s) / (1 - eps_MC) = a + eps_d * (1 - s) / (1 - eps_MC)
//
// so the extra factor on top of the SF already applied is
//
//   a'/a = 1 + eps_d * (1 - s) / (1 - eps_d),    1 - eps_d = a * (1 - eps_MC).
//
// Written this way it needs only eps_MC^veto -- the MC-truth veto efficiency
// map the anti-veto SF was itself built from, see smoothVetoSF.py -- plus the
// SF value being applied, so it stays consistent with whichever veto SF
// flavour the analysis runs with (--useRefinedVeto or not), and reduces to
// (1 - eps_MC * s)/(1 - eps_MC) when no veto SF is applied at all (a = 1).
//
// A word of warning on the size. 1 - eps_d is only ~1% above pt ~ 30 GeV, so
// even a few percent of refit inefficiency dominates the probability of
// escaping the veto and the factor reaches O(10) over much of the hotspots,
// O(100) in the worst cell. That is not a defect of the formula: in those
// cells the refit hole really is the main way a second muon goes unvetoed in
// data. It does mean the migrated background is modelled by up-weighting the
// few MC events that fail the veto for ordinary reasons, which is where this
// approach pays a statistical price relative to simulating the migration.
//
// Two known limitations, both confined to the hotspot cells:
//
//  * the veto SF stat/syst variations are applied as a ratio to the whole
//    nominal weight, i.e. as a * f -> a_var * f, whereas the correct varied
//    weight is a_var + eps_d*(1-s)/(1-eps_MC), the second term being unchanged.
//    The veto SF nuisances are therefore inflated by up to the size of the
//    factor where it is large. Conservative, but not a faithful variation.
//
//  * no uncertainty is assigned to the refit SF map itself, here or in
//    cvh_muon_sf.
template <typename HIST_EFF> class cvh_antiveto_helper {

public:
  cvh_antiveto_helper(HIST_EFF &&eff_veto_mc)
      : eff_(std::make_shared<const HIST_EFF>(std::move(eff_veto_mc))) {}

  // pt, eta, phi, charge of the postFSR gen muon that failed the veto, with
  // charge = -99 flagging 'no such muon in acceptance' as for the veto SF
  // helpers; antivetoSF is the anti-veto SF already multiplied into the event
  // weight for that same muon (pass 1 if none is).
  double operator()(float pt, float eta, float phi, int charge,
                    double antivetoSF) const {

    if (charge <= -99)
      return 1.0;

    // gen quantities stand in for the track that was never reconstructed; the
    // charge and pt only serve to undo the bending to the module radius
    const double sf = cvh_muon_sf(eta, phi, charge, pt);
    if (sf >= 1.0)
      return 1.0; // outside the hotspots, i.e. almost always

    const double failMC = 1.0 - veto_efficiency_mc(pt, eta, charge);
    if (!(failMC > 0.0))
      return 1.0; // empty or pathological bin of the efficiency map

    const double failData = antivetoSF * failMC; // = 1 - eps_d
    if (!(failData > 0.0) || failData >= 1.0)
      return 1.0;

    return 1.0 + (1.0 - failData) * (1.0 - sf) / failData;
  }

private:
  double veto_efficiency_mc(float pt, float eta, int charge) const {
    // the map has no flow bins, and the gen muon can sit above its pt range
    const int ieta = clamped_index<0>(eta);
    const int ipt = clamped_index<1>(pt);
    const int icharge = clamped_index<2>(charge);
    return eff_->at(ieta, ipt, icharge).value();
  }

  template <unsigned int Iaxis, typename T> int clamped_index(T value) const {
    const auto &axis = eff_->template axis<Iaxis>();
    return std::clamp(axis.index(value), 0, static_cast<int>(axis.size()) - 1);
  }

  std::shared_ptr<const HIST_EFF> eff_;
};

} // namespace wrem

#endif
