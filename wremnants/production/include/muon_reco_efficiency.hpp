#ifndef WREMNANTS_MUON_RECO_EFFICIENCIES_H
#define WREMNANTS_MUON_RECO_EFFICIENCIES_H

#include <TH3D.h>
#include <algorithm>
#include <memory>
#include <vector>

#include "utils.hpp"

namespace wrem {

class muon_reco_efficiency_helper {
public:
  muon_reco_efficiency_helper(const TH3D &input_hist)
      : sf_hist(std::make_shared<TH3D>(input_hist)) {}

  double get_weight(float eta, float pt) const {
    int binX = sf_hist->GetXaxis()->FindBin(eta);
    int binY = sf_hist->GetYaxis()->FindBin(pt);

    binX = std::max(1, std::min(binX, sf_hist->GetNbinsX()));
    binY = std::max(1, std::min(binY, sf_hist->GetNbinsY()));

    return sf_hist->GetBinContent(binX, binY, 1);
  }

  // Vector version --> i dont know why this is needed because the other helpers
  // seem to just have the float version
  double operator()(const Vec_f &etas, const Vec_f &pts) const {
    double weight = 1.0;
    for (size_t i = 0; i < etas.size(); ++i) {
      weight *= get_weight(etas[i], pts[i]);
    }
    return weight;
  }

private:
  std::shared_ptr<TH3D> sf_hist;
};

} // namespace wrem

#endif