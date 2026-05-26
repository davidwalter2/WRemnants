#ifndef WREMNANTS_MUON_EFFICIENCIES_SMOOTH_H
#define WREMNANTS_MUON_EFFICIENCIES_SMOOTH_H

#include "defines.hpp"
#include <array>
#include <boost/histogram/axis.hpp>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

namespace wrem {

// TODO use enums for integer/boolean/category axes so that the code is less
// error-prone?

template <int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_base {
public:
  muon_efficiency_smooth_helper_base(HIST_SF &&sf_all)
      : sf_all_(std::make_shared<const HIST_SF>(std::move(sf_all))) {}

  std::array<double, 2> scale_factor_array(int pt_idx, int eta_idx,
                                           int sapt_idx, int saeta_idx,
                                           int charge_idx,
                                           int idx_nom_alt) const {

    auto const eff_type_idx_reco = idx_reco_;
    auto const eff_type_idx_tracking = idx_tracking_;

    const double reco =
        sf_all_->at(eta_idx, pt_idx, charge_idx, eff_type_idx_reco, idx_nom_alt)
            .value();
    const double tracking = sf_all_
                                ->at(saeta_idx, sapt_idx, charge_idx,
                                     eff_type_idx_tracking, idx_nom_alt)
                                .value();

    std::array<double, 2> ret = {reco, tracking};

    return ret;
  }

  double scale_factor_product(float pt, float eta, float sapt, float saeta,
                              int charge, int idx_nom_alt) const {

    auto const eta_idx = sf_all_->template axis<0>().index(eta);
    auto const pt_idx = sf_all_->template axis<1>().index(pt);
    auto const charge_idx = sf_all_->template axis<2>().index(charge);
    auto const saeta_idx = sf_all_->template axis<0>().index(saeta);
    auto const sapt_idx = sf_all_->template axis<1>().index(sapt);

    std::array<double, 2> allSF = scale_factor_array(
        pt_idx, eta_idx, sapt_idx, saeta_idx, charge_idx, idx_nom_alt);
    double sf = 1.0;
    for (int i = 0; i < allSF.size(); i++) {
      sf *= allSF[i];
    }
    return sf;
  }

  using syst_tensor_t =
      Eigen::TensorFixedSize<double, Eigen::Sizes<2, NSysts>>; // 2 bins for
                                                               // reco, tracking

  syst_tensor_t sf_syst_var(float pt, float eta, float sapt, float saeta,
                            int charge) const {

    syst_tensor_t res;

    auto const eta_idx = sf_all_->template axis<0>().index(eta);
    auto const pt_idx = sf_all_->template axis<1>().index(pt);
    auto const saeta_idx = sf_all_->template axis<0>().index(saeta);
    auto const sapt_idx = sf_all_->template axis<1>().index(sapt);
    auto const charge_idx = sf_all_->template axis<2>().index(charge);

    std::array<double, 2> allSF_nomi = scale_factor_array(
        pt_idx, eta_idx, sapt_idx, saeta_idx, charge_idx, idx_nom_);

    for (int ns = 0; ns < NSysts; ns++) {

      std::array<double, 2> allSF_alt = scale_factor_array(
          pt_idx, eta_idx, sapt_idx, saeta_idx, charge_idx,
          sf_all_->template axis<4>().index(
              ns + 1)); // 0 is the nominal, systs starts from 1

      // anticorrelation between iso and antiiso already embedded in the numbers
      // stored in the histograms also the alternate comes from data efficiency
      // variation only, so the anticorrelation in the efficiencies is preserved
      // in the scale factors

      // order is reco-tracking-idip-trigger-iso
      for (int i = 0; i < allSF_nomi.size(); i++) {
        // if (allSF_nomi[i] <= 0)
        //     std::cout << "allSF_nomi/alt[" << i << "] = " << allSF_nomi[i] <<
        //     "/" << allSF_alt[i] << " --> pt/eta/charge = " << pt << "/" <<
        //     eta << "/" << charge << std::endl;
        res(i, ns) = allSF_alt[i] / allSF_nomi[i];
      }
    }

    return res;
  }

protected:
  std::shared_ptr<const HIST_SF> sf_all_;
  // cache the bin indices since the string category lookup is slow
  int idx_reco_ = sf_all_->template axis<3>().index("reco");
  int idx_tracking_ = sf_all_->template axis<3>().index("tracking");

  int idx_nom_ = sf_all_->template axis<4>().index(0);
};

// base template for one-lepton case
template <AnalysisType analysisType, int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper
    : public muon_efficiency_smooth_helper_base<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base<NSysts, HIST_SF>;
  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper(const base_t &other) : base_t(other) {}

  double operator()(float pt, float eta, float sapt, float saeta, int charge) {
    return base_t::scale_factor_product(pt, eta, sapt, saeta, charge,
                                        base_t::idx_nom_);
  }
};

// specialization for two-lepton case Wlike
template <int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper<AnalysisType::Wlike, NSysts, HIST_SF>
    : public muon_efficiency_smooth_helper_base<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base<NSysts, HIST_SF>;
  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper(const base_t &other) : base_t(other) {}

  double operator()(float trig_pt, float trig_eta, float trig_sapt,
                    float trig_saeta, int trig_charge, float nontrig_pt,
                    float nontrig_eta, float nontrig_sapt, float nontrig_saeta,
                    int nontrig_charge) {
    const double sftrig =
        base_t::scale_factor_product(trig_pt, trig_eta, trig_sapt, trig_saeta,
                                     trig_charge, base_t::idx_nom_);
    const double sfnontrig = base_t::scale_factor_product(
        nontrig_pt, nontrig_eta, nontrig_sapt, nontrig_saeta, nontrig_charge,
        base_t::idx_nom_);
    return sftrig * sfnontrig;
  }
};

// specialization for two-lepton case Dilepton
template <int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper<AnalysisType::Dilepton, NSysts, HIST_SF>
    : public muon_efficiency_smooth_helper_base<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base<NSysts, HIST_SF>;
  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper(const base_t &other) : base_t(other) {}

  // may also assume the first is the one passing the trigger for sure, but it
  // depends on how these operators are called in the loop keeping both flags is
  // redundant but more flexible since there is no assumption on the sorting
  double operator()(float first_pt, float first_eta, float first_sapt,
                    float first_saeta, int first_charge, float second_pt,
                    float second_eta, float second_sapt, float second_saeta,
                    int second_charge) {

    const double sftrig = base_t::scale_factor_product(
        first_pt, first_eta, first_sapt, first_saeta, first_charge,
        base_t::idx_nom_);
    const double sfnontrig = base_t::scale_factor_product(
        second_pt, second_eta, second_sapt, second_saeta, second_charge,
        base_t::idx_nom_);
    return sftrig * sfnontrig;
  }
};

// Now the syst, which is similar to the nominal
//
// base template for one lepton case
template <AnalysisType analysisType, int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_syst
    : public muon_efficiency_smooth_helper_base<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base<NSysts, HIST_SF>;
  using tensor_t = typename base_t::syst_tensor_t;

  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper_syst(const base_t &other) : base_t(other) {}

  tensor_t operator()(float pt, float eta, float sapt, float saeta, int charge,
                      bool pass_iso, double nominal_weight = 1.0) {
    constexpr bool iso_with_trigger = true;
    constexpr bool pass_trigger = true;
    return nominal_weight * base_t::sf_syst_var(pt, eta, sapt, saeta, charge,
                                                pass_iso, pass_trigger,
                                                iso_with_trigger);
  }
};

// specialization for two-lepton case Wlike
template <int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_syst<AnalysisType::Wlike, NSysts, HIST_SF>
    : public muon_efficiency_smooth_helper_base<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base<NSysts, HIST_SF>;
  using tensor_t = typename base_t::syst_tensor_t;

  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper_syst(const base_t &other) : base_t(other) {}

  tensor_t operator()(float trig_pt, float trig_eta, float trig_sapt,
                      float trig_saeta, int trig_charge, float nontrig_pt,
                      float nontrig_eta, float nontrig_sapt,
                      float nontrig_saeta, int nontrig_charge,
                      double nominal_weight = 1.0) {
    const tensor_t variation_trig = base_t::sf_syst_var(
        trig_pt, trig_eta, trig_sapt, trig_saeta, trig_charge);
    const tensor_t variation_nontrig = base_t::sf_syst_var(
        nontrig_pt, nontrig_eta, nontrig_sapt, nontrig_saeta, nontrig_charge);
    return nominal_weight * variation_trig * variation_nontrig;
  }
};

// specialization for two-lepton case Dilepton
template <int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_syst<AnalysisType::Dilepton, NSysts,
                                         HIST_SF>
    : public muon_efficiency_smooth_helper_base<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base<NSysts, HIST_SF>;
  using tensor_t = typename base_t::syst_tensor_t;

  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper_syst(const base_t &other) : base_t(other) {}

  // may also assume the first is the one passing the trigger for sure, but it
  // depends on how these operators are called in the loop keeping both flags is
  // redundant but more flexible since there is no assumption on the sorting
  tensor_t operator()(float first_pt, float first_eta, float first_sapt,
                      float first_saeta, int first_charge, float second_pt,
                      float second_eta, float second_sapt, float second_saeta,
                      int second_charge, double nominal_weight = 1.0) {
    const tensor_t variation_trig = base_t::sf_syst_var(
        first_pt, first_eta, first_sapt, first_saeta, first_charge);
    const tensor_t variation_nontrig = base_t::sf_syst_var(
        second_pt, second_eta, second_sapt, second_saeta, second_charge);
    return nominal_weight * variation_trig * variation_nontrig;
  }
};

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

////
// STAT UNCERTAINTY
// this is now an independent class with respect to the previous one which only
// deals with nominal and statistical variations
////
//// BASE CLASS FOR HELPER_STAT
template <int NEtaBins, int NPtEigenBins, int NCharges, typename HIST_SF>
class muon_efficiency_smooth_helper_stat_base {
public:
  muon_efficiency_smooth_helper_stat_base(HIST_SF &&sf_type)
      : sf_type_(std::make_shared<const HIST_SF>(std::move(sf_type))) {}

  // number of eta bins, number of eigen variations for pt axis, then 2 charges
  using stat_tensor_t =
      Eigen::TensorFixedSize<double,
                             Eigen::Sizes<NEtaBins, NPtEigenBins, NCharges>>;

  int checkEffTypeInAxis(boost::histogram::axis::category<std::string> axis,
                         const std::string &match = "match") {
    int ret = -1;
    for (Int_t i = 0; i < axis.size(); i++) {
      if (match == axis.value(i)) {
        ret = i;
        break;
      }
    }
    return ret;
  }

  // general case with no isolation (and no uT dependence)
  stat_tensor_t sf_stat_var(float pt, float eta, int charge) const {
    stat_tensor_t res;
    res.setConstant(1.0);

    auto const eta_idx = sf_type_->template axis<0>().index(eta);
    auto const pt_idx = sf_type_->template axis<1>().index(pt);
    auto const charge_idx = sf_type_->template axis<2>().index(charge);
    auto const eff_type_idx = 0;
    auto const eigen_axis = sf_type_->template axis<4>();

    // overflow/underflow are attributed to adjacent bin
    auto const tensor_eta_idx = std::clamp(eta_idx, 0, NEtaBins - 1);

    auto const &cell_nomi =
        sf_type_->at(eta_idx, pt_idx, charge_idx, eff_type_idx, idx_nom_);

    const double sf_nomi = cell_nomi.value();
    // loop on dimension with the parameters variations (the histogram already
    // has the alternate SF) start from 1 because first bin contains the nominal
    // SF
    for (int tensor_eigen_idx = 1; tensor_eigen_idx <= NPtEigenBins;
         tensor_eigen_idx++) {

      auto const eigen_axis_idx = eigen_axis.index(tensor_eigen_idx);

      auto const &cell_stat = sf_type_->at(eta_idx, pt_idx, charge_idx,
                                           eff_type_idx, eigen_axis_idx);

      const double sf_stat = cell_stat.value();
      const double sf_stat_variation = sf_stat / sf_nomi;

      res(tensor_eta_idx, tensor_eigen_idx - 1, charge_idx) *=
          sf_stat_variation;
    }

    return res;
  }

  // special case for isolation (and no uT dependence)
  stat_tensor_t sf_stat_var_iso(float pt, float eta, int charge) const {
    stat_tensor_t res;
    res.setConstant(1.0);

    auto const eta_idx = sf_type_->template axis<0>().index(eta);
    auto const pt_idx = sf_type_->template axis<1>().index(pt);
    auto const charge_idx = sf_type_->template axis<2>().index(charge);
    auto const eigen_axis = sf_type_->template axis<4>();

    // overflow/underflow are attributed to adjacent bin
    auto const tensor_eta_idx = std::clamp(eta_idx, 0, NEtaBins - 1);

    auto const &cell_nomi = sf_type_->at(eta_idx, pt_idx, charge_idx, idx_nom_);

    const double sf_nomi = cell_nomi.value();
    // loop on dimension with the parameters variations (the histogram already
    // has the alternate SF) start from 1 because first bin contains the nominal
    // SF
    for (int tensor_eigen_idx = 1; tensor_eigen_idx <= NPtEigenBins;
         tensor_eigen_idx++) {

      auto const eigen_axis_idx = eigen_axis.index(tensor_eigen_idx);

      auto const &cell_stat =
          sf_type_->at(eta_idx, pt_idx, charge_idx, eigen_axis_idx);

      const double sf_stat = cell_stat.value();
      const double sf_stat_variation = sf_stat / sf_nomi;
      res(tensor_eta_idx, tensor_eigen_idx - 1, charge_idx) *=
          sf_stat_variation;
    }

    return res;
  }

  // special case for isolation with ut dependence
  stat_tensor_t sf_stat_var_iso(float pt, float eta, float ut,
                                int charge) const {
    stat_tensor_t res;
    res.setConstant(1.0);

    auto const eta_idx = sf_type_->template axis<0>().index(eta);
    auto const pt_idx = sf_type_->template axis<1>().index(pt);
    auto const charge_idx = sf_type_->template axis<2>().index(charge);
    auto const eigen_axis = sf_type_->template axis<4>();

    // overflow/underflow are attributed to adjacent bin
    auto const tensor_eta_idx = std::clamp(eta_idx, 0, NEtaBins - 1);

    auto const &cell_nomi = sf_type_->at(eta_idx, pt_idx, charge_idx, idx_nom_);

    const double sf_nomi = cell_nomi.value();
    // loop on dimension with the parameters variations (the histogram already
    // has the alternate SF) start from 1 because first bin contains the nominal
    // SF
    for (int tensor_eigen_idx = 1; tensor_eigen_idx <= NPtEigenBins;
         tensor_eigen_idx++) {

      auto const eigen_axis_idx = eigen_axis.index(tensor_eigen_idx);

      auto const &cell_stat =
          sf_type_->at(eta_idx, pt_idx, charge_idx, eigen_axis_idx);

      const double sf_stat = cell_stat.value();
      const double sf_stat_variation = sf_stat / sf_nomi;
      res(tensor_eta_idx, tensor_eigen_idx - 1, charge_idx) *=
          sf_stat_variation;
    }

    return res;
  }

protected:
  std::shared_ptr<const HIST_SF> sf_type_;
  // cache the bin indices since the string category lookup is slow
  int idx_nom_ = sf_type_->template axis<4>().index(0);
  // input effStat axis is organized as nomi - UpVar, with nomi centered at 0
};

////
//// General case no uT dependence (isolation is treated separately)
////

// base template for one lepton case
template <AnalysisType analysisType, int NEtaBins, int NPtEigenBins,
          int NCharges, typename HIST_SF>
class muon_efficiency_smooth_helper_stat
    : public muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins,
                                                     NCharges, HIST_SF> {

public:
  using stat_base_t =
      muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges,
                                              HIST_SF>;
  using tensor_t = typename stat_base_t::stat_tensor_t;

  using stat_base_t::stat_base_t;

  tensor_t operator()(float pt, float eta, int charge,
                      double nominal_weight = 1.0) {
    return nominal_weight * stat_base_t::sf_stat_var(pt, eta, charge);
  }
};

// specialization for two-lepton case Dilepton
template <int NEtaBins, int NPtEigenBins, int NCharges, typename HIST_SF>
class muon_efficiency_smooth_helper_stat<AnalysisType::Dilepton, NEtaBins,
                                         NPtEigenBins, NCharges, HIST_SF>
    : public muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins,
                                                     NCharges, HIST_SF> {

public:
  using stat_base_t =
      muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges,
                                              HIST_SF>;
  using tensor_t = typename stat_base_t::stat_tensor_t;

  using stat_base_t::stat_base_t;

  tensor_t operator()(float first_pt, float first_eta, int first_charge,
                      float second_pt, float second_eta, int second_charge,
                      double nominal_weight = 1.0) {

    const tensor_t variation_first =
        stat_base_t::sf_stat_var(first_pt, first_eta, first_charge);
    const tensor_t variation_second =
        stat_base_t::sf_stat_var(second_pt, second_eta, second_charge);
    return nominal_weight * variation_first * variation_second;
  }
};

// base template for one lepton case
template <AnalysisType analysisType, int NEtaBins, int NPtEigenBins,
          int NCharges, typename HIST_SF>
class muon_efficiency_smooth_helper_stat_utDep
    : public muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins,
                                                     NCharges, HIST_SF> {

public:
  using stat_base_t =
      muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges,
                                              HIST_SF>;
  using tensor_t = typename stat_base_t::stat_tensor_t;

  using stat_base_t::stat_base_t;

  tensor_t operator()(float pt, float eta, float ut, int charge,
                      double nominal_weight = 1.0) {
    constexpr bool pass_trigger = true;
    return nominal_weight *
           stat_base_t::sf_stat_var(pt, eta, ut, charge, pass_trigger);
  }
};

// specialization for two-lepton case Dilepton
template <int NEtaBins, int NPtEigenBins, int NCharges, typename HIST_SF>
class muon_efficiency_smooth_helper_stat_utDep<AnalysisType::Dilepton, NEtaBins,
                                               NPtEigenBins, NCharges, HIST_SF>
    : public muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins,
                                                     NCharges, HIST_SF> {

public:
  using stat_base_t =
      muon_efficiency_smooth_helper_stat_base<NEtaBins, NPtEigenBins, NCharges,
                                              HIST_SF>;
  using tensor_t = typename stat_base_t::stat_tensor_t;

  using stat_base_t::stat_base_t;

  tensor_t operator()(float first_pt, float first_eta, int first_charge,
                      float second_pt, float second_eta, int second_charge,
                      double nominal_weight = 1.0) {

    const tensor_t variation_first =
        stat_base_t::sf_stat_var(first_pt, first_eta, first_charge);
    const tensor_t variation_second =
        stat_base_t::sf_stat_var(second_pt, second_eta, second_charge);
    return nominal_weight * variation_first * variation_second;
  }
};

///////////
// ADDITIONAL UTILITY HELPERS TO MANAGE A SINGLE EFFICIENCY STEP
// originally conceived to deal with additional systematic uncertainties on the
// tracking step the "nominal" helper is also implemented so one can add special
// additional scale factors if ever useful
///////////
template <int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_base_oneStep {
public:
  muon_efficiency_smooth_helper_base_oneStep(HIST_SF &&sf_type)
      : sf_type_(std::make_shared<const HIST_SF>(std::move(sf_type))) {}

  double scale_factor_byCell(int pt_idx, int eta_idx, int charge_idx,
                             int idx_nom_alt) const {
    const double sf = sf_type_->at(eta_idx, pt_idx, charge_idx, idx_nom_alt);
    // std::cout << "Scale factor " << sf << std::endl;
    return sf;
  }

  double scale_factor(float pt, float eta, int charge, int idx_nom_alt) const {

    auto const eta_idx = sf_type_->template axis<0>().index(eta);
    auto const pt_idx = sf_type_->template axis<1>().index(pt);
    auto const charge_idx = sf_type_->template axis<2>().index(charge);
    return scale_factor_byCell(pt_idx, eta_idx, charge_idx, idx_nom_alt);
  }

  using syst_tensor_t = Eigen::TensorFixedSize<double, Eigen::Sizes<NSysts>>;

  syst_tensor_t sf_syst_var(float pt, float eta, int charge) const {

    syst_tensor_t res;

    auto const eta_idx = sf_type_->template axis<0>().index(eta);
    auto const pt_idx = sf_type_->template axis<1>().index(pt);
    auto const charge_idx = sf_type_->template axis<2>().index(charge);

    double sf_nomi = scale_factor_byCell(pt_idx, eta_idx, charge_idx, idx_nom_);

    for (int ns = 0; ns < NSysts; ns++) {

      double sf_alt = scale_factor_byCell(
          pt_idx, eta_idx, charge_idx,
          sf_type_->template axis<3>().index(
              ns + 1)); // 0 is the nominal, systs starts from 1
      res(ns) = sf_alt / sf_nomi;
    }

    return res;
  }

protected:
  std::shared_ptr<const HIST_SF> sf_type_;
  // cache the bin indices since the string category lookup is slow
  int idx_nom_ = sf_type_->template axis<3>().index(
      0); // input effStat axis is organized as nomi - UpVar, with nomi centered
          // at 0
};

// base template for one-lepton case
template <AnalysisType analysisType, int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_oneStep
    : public muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF>;
  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper_oneStep(const base_t &other) : base_t(other) {}

  double operator()(float pt, float eta, int charge) {
    return base_t::scale_factor(pt, eta, charge, base_t::idx_nom_);
  }
};

// specialization for two-lepton case Wlike
template <int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_oneStep<AnalysisType::Wlike, NSysts,
                                            HIST_SF>
    : public muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF>;
  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper_oneStep(const base_t &other) : base_t(other) {}

  double operator()(float trig_pt, float trig_eta, int trig_charge,
                    float nontrig_pt, float nontrig_eta, int nontrig_charge) {
    const double sftrig =
        base_t::scale_factor(trig_pt, trig_eta, trig_charge, base_t::idx_nom_);
    const double sfnontrig = base_t::scale_factor(
        nontrig_pt, nontrig_eta, nontrig_charge, base_t::idx_nom_);
    return sftrig * sfnontrig;
  }
};

// specialization for two-lepton case Dilepton
template <int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_oneStep<AnalysisType::Dilepton, NSysts,
                                            HIST_SF>
    : public muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF>;
  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper_oneStep(const base_t &other) : base_t(other) {}

  // may also assume the first is the one passing the trigger for sure, but it
  // depends on how these operators are called in the loop keeping both flags is
  // redundant but more flexible since there is no assumption on the sorting
  double operator()(float first_pt, float first_eta, int first_charge,
                    float second_pt, float second_eta, int second_charge) {
    const double sftrig = base_t::scale_factor(first_pt, first_eta,
                                               first_charge, base_t::idx_nom_);
    const double sfnontrig = base_t::scale_factor(
        second_pt, second_eta, second_charge, base_t::idx_nom_);
    return sftrig * sfnontrig;
  }
};

// Now the syst, which is similar to the nominal
//
// base template for one lepton case
template <AnalysisType analysisType, int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_syst_oneStep
    : public muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF>;
  using tensor_t = typename base_t::syst_tensor_t;

  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper_syst_oneStep(const base_t &other)
      : base_t(other) {}

  tensor_t operator()(float pt, float eta, int charge,
                      double nominal_weight = 1.0) {
    return nominal_weight * base_t::sf_syst_var(pt, eta, charge);
  }
};

// specialization for two-lepton case Wlike
template <int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_syst_oneStep<AnalysisType::Wlike, NSysts,
                                                 HIST_SF>
    : public muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF>;
  using tensor_t = typename base_t::syst_tensor_t;

  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper_syst_oneStep(const base_t &other)
      : base_t(other) {}

  tensor_t operator()(float trig_pt, float trig_eta, int trig_charge,
                      float nontrig_pt, float nontrig_eta, int nontrig_charge,
                      double nominal_weight = 1.0) {
    const tensor_t variation_trig =
        base_t::sf_syst_var(trig_pt, trig_eta, trig_charge);
    const tensor_t variation_nontrig =
        base_t::sf_syst_var(nontrig_pt, nontrig_eta, nontrig_charge);
    return nominal_weight * variation_trig * variation_nontrig;
  }
};

// specialization for two-lepton case Dilepton
template <int NSysts, typename HIST_SF>
class muon_efficiency_smooth_helper_syst_oneStep<AnalysisType::Dilepton, NSysts,
                                                 HIST_SF>
    : public muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF> {

public:
  using base_t = muon_efficiency_smooth_helper_base_oneStep<NSysts, HIST_SF>;
  using tensor_t = typename base_t::syst_tensor_t;

  // inherit constructor
  using base_t::base_t;

  muon_efficiency_smooth_helper_syst_oneStep(const base_t &other)
      : base_t(other) {}

  // may also assume the first is the one passing the trigger for sure, but it
  // depends on how these operators are called in the loop keeping both flags is
  // redundant but more flexible since there is no assumption on the sorting
  tensor_t operator()(float first_pt, float first_eta, int first_charge,
                      float second_pt, float second_eta, int second_charge,
                      double nominal_weight = 1.0) {
    constexpr bool iso_with_trigger =
        true; // will be P(iso|passTrigger) or P(iso|failTrigger) depending on
              // first_passtrigger and second_passtrigger
    constexpr bool pass_iso = true;
    const tensor_t variation_trig =
        base_t::sf_syst_var(first_pt, first_eta, first_charge);
    const tensor_t variation_nontrig =
        base_t::sf_syst_var(second_pt, second_eta, second_charge);
    return nominal_weight * variation_trig * variation_nontrig;
  }
};

} // namespace wrem

#endif
