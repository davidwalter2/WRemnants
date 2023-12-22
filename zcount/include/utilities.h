#ifndef WREMNANTS_UTILITIES_H
#define WREMNANTS_UTILITIES_H

#include "definitions.h"

using namespace ROOT;

namespace zcount {

double deltaPhi(float phi1, float phi2) {
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2.0*M_PI;
    while (result <= -1.0*M_PI) result += 2.0*M_PI;
    return result;
}

double deltaR2(float eta1, float phi1, float eta2, float phi2) {
    double deta = eta1-eta2;
    double dphi = deltaPhi(phi1,phi2);
    return deta*deta + dphi*dphi;
}

Vec_b hasMatchDR2(const Vec_f& vec1_eta, const Vec_f& vec1_phi, const float eta, const float phi, const float dr2) {
    Vec_b res(vec1_eta.size(), false); // initialize to 0
    for (unsigned int i = 0; i < res.size(); ++i) {
        if (deltaR2(vec1_eta[i], vec1_phi[i], eta, phi) < dr2){
            res[i] = true;
        }
    }
    return res;
}

Vec_b hasMatchDR2(const Vec_f& vec1_eta, const Vec_f& vec1_phi, const Vec_f& vec2_eta, const Vec_f& vec2_phi, const float dr2) {
    // make one to one comparisons
    Vec_b res(vec1_eta.size(), false); // initialize to 0
    for (unsigned int i = 0; i < res.size(); ++i) {
        if (deltaR2(vec1_eta[i], vec1_phi[i], vec2_eta[i], vec2_phi[i]) < dr2){
            res[i] = true;
        }
    }
    return res;
}

bool hasAnyMatchDR2(const float eta, const float phi, const Vec_f& vec_eta, const Vec_f& vec_phi, const float dr2) {
    for (unsigned int i = 0; i < vec_eta.size(); ++i) {
        if (deltaR2(eta, phi, vec_eta[i], vec_phi[i]) < dr2){
            return true;
        }
    }
    return false;
}

Vec_b hasAnyMatchDR2(const Vec_f& vec1_eta, const Vec_f& vec1_phi, const Vec_f& vec2_eta, const Vec_f& vec2_phi, const float dr2) {
    Vec_b res(vec1_eta.size(), false); // initialize to 0
    for (unsigned int i = 0; i < res.size(); ++i) {
        for (unsigned int j = 0; j < vec2_eta.size(); ++j) {
            if (deltaR2(vec1_eta[i], vec1_phi[i], vec2_eta[j], vec2_phi[j]) < dr2){
                res[i] = true;
                break;
            }
        }
    }
    return res;
}

Vec_lorentz makeLorentzVector(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const float mass) {
    // std::cout<<"Make lorentzvector"<<std::endl;
    Vec_lorentz res(pt.size());
    for (unsigned int i = 0; i < res.size(); ++i) {
        res[i] = ROOT::Math::PtEtaPhiMVector(pt[i], eta[i], phi[i], mass);
        // std::cout<<"--- "<< i <<std::endl;
        // std::cout<<"pT = "<<pt[i] <<std::endl;
        // std::cout<<"eta = "<<eta[i] <<std::endl;
        // std::cout<<"phi = "<<phi[i] <<std::endl;
    }
    // std::cout<<"End Make lorentzvector"<<std::endl;
    return res;
}

// --- for systematics
Eigen::TensorFixedSize<double, Eigen::Sizes<2>> constantScaling(double nominal_weight, double scale) {
    Eigen::TensorFixedSize<double, Eigen::Sizes<2>> outWeights;
    // Down weight, then up weight
    outWeights(0) = nominal_weight / scale;
    outWeights(1) = nominal_weight * scale;
    return outWeights;
}

Eigen::TensorFixedSize<double, Eigen::Sizes<2>> twoPointScaling(double nominal_weight, double scaleDown, double scaleUp) {
    Eigen::TensorFixedSize<double, Eigen::Sizes<2>> outWeights;
    // Down weight, then up weight, nominal_weight should not already include a centralScale weight if any
    outWeights(0) = nominal_weight * scaleDown;
    outWeights(1) = nominal_weight * scaleUp;
    return outWeights;
}

// --- for selection
Vec_b mergedStandAloneIsGoodGlobal(
    const Vec_i& MergedStandAloneMuon_extraIdx, 
    const Vec_f& MergedStandAloneMuon_eta, 
    const Vec_f& MergedStandAloneMuon_phi, 
    const Vec_i& Muon_standaloneExtraIdx, 
    const Vec_f& Muon_eta, 
    const Vec_f& Muon_phi
){
    Vec_b res(MergedStandAloneMuon_extraIdx.size(), false);
    for (unsigned int i = 0; i < res.size(); ++i) {
        for (unsigned int j = 0; j < Muon_standaloneExtraIdx.size(); ++j) {
            if(MergedStandAloneMuon_extraIdx[i] == Muon_standaloneExtraIdx[j]){
                // found corresponding muon object
                if(deltaR2(MergedStandAloneMuon_eta[i], MergedStandAloneMuon_phi[i], Muon_eta[j], Muon_phi[j]) < 0.09
                ){
                    res[i] = true;
                }
                break;
            }
        }
    }
    return res;
}

Vec_b hasTrackRef(
    const Vec_i& extraIdx, 
    const Vec_i& reference_extraIdx
){
    Vec_b res(extraIdx.size(), false);
    for (unsigned int i = 0; i < res.size(); ++i) {
        for (unsigned int j = 0; j < reference_extraIdx.size(); ++j) {
            if(extraIdx[i] == reference_extraIdx[j]){
                // found corresponding reference object
                res[i] = true;
                break;
            }
        }
    }
    return res;    
}

}


#endif