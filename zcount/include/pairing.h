#ifndef WREMNANTS_PAIRING_H
#define WREMNANTS_PAIRING_H

#include "definitions.h"
#include "utilities.h"

using namespace ROOT;

namespace zcount {

using LorentzVector = ROOT::Math::PtEtaPhiMVector;

using Vec_lorentz = ROOT::VecOps::RVec<LorentzVector>;

// --- pair muons and calculate mass
Vec_f makePairsMass(const Vec_lorentz& l1) {
    if(l1.size() <= 1) {
        Vec_f res(0);
        return res;
    }
    const unsigned int combinations = l1.size()*(l1.size()-1)/2;
    Vec_f res(combinations); 
    for (unsigned int i = 0, k=0; i < l1.size(); ++i) {
        for (unsigned int j = i+1; j < l1.size(); ++j) {
            res[k] = (l1[i] + l1[j]).mass();
            k++;
        }
    }
    return res;
}
Vec_f makePairsMass(const Vec_lorentz& l1, const Vec_lorentz& l2) {
    // std::cout<<"Make pairs mass"<<std::endl;
    const unsigned int combinations = l1.size()*l2.size();
    Vec_f res(combinations);
    for (unsigned int i = 0; i < l1.size(); ++i) {
        for (unsigned int j = 0; j < l2.size(); ++j) {
            res[i*l2.size()+j] = (l1[i] + l2[j]).mass();
        }
    }
    // std::cout<<"End Make pairs mass"<<std::endl;
    return res;
}

Vec_f makePairsMass(const Vec_lorentz& l1, const LorentzVector l2) {
    Vec_f res(l1.size());
    for (unsigned int i = 0; i < l1.size(); ++i) {
        res[i] = (l1[i] + l2).mass();
    }
    return res;
}


double makePairsMass(const LorentzVector l1, const LorentzVector l2) {
    return (l1 + l2).mass();
}

// --- pair muons and check if two muons have opposite charge
Vec_b makePairsOS(const Vec_i& l1) {
    if(l1.size() <= 1) {
        Vec_b res(0);
        return res;
    }
    const unsigned int combinations = l1.size()*(l1.size()-1)/2;
    Vec_b res(combinations);
    for (unsigned int i = 0, k=0; i < l1.size(); ++i) {
        for (unsigned int j = i+1; j < l1.size(); ++j) {
            res[k] = (l1[i] != l1[j]);
            k++;
        }
    }
    return res;
}

Vec_b makePairsOS(const Vec_i& l1, const Vec_i& l2) {
    // std::cout<<"Make pairs OS"<<std::endl;
    const unsigned int combinations = l1.size()*l2.size();
    Vec_b res(combinations);
    for (unsigned int i = 0; i < l1.size(); ++i) {
        for (unsigned int j = 0; j < l2.size(); ++j) {
            res[i*l2.size()+j] = (l1[i] != l2[j]);
        }
    }
    // std::cout<<"End Make pairs OS"<<std::endl;
    return res;
}

Vec_b makePairsOS(const Vec_i& l1, const int l2) {
    Vec_b res(l1.size());
    for (unsigned int i = 0; i < l1.size(); ++i) {
        res[i] = l1[i] != l2;
    }
    return res;
}

bool makePairsOS(const int c1, const int c2) {
    return c1 != c2;
}

// --- pair muons and check if two muons are both gen matched 
Vec_b makePairsGen(const Vec_b& l1) {
    if(l1.size() <= 1) {
        Vec_b res(0);
        return res;
    }
    const unsigned int combinations = l1.size()*(l1.size()-1)/2;
    Vec_b res(combinations);
    for (unsigned int i = 0, k=0; i < l1.size(); ++i) {
        for (unsigned int j = i+1; j < l1.size(); ++j) {
            res[k] = l1[i] & l1[j];
            k++;
        }
    }
    return res;
}

Vec_b makePairsGen(const Vec_b& l1, const Vec_b& l2) {
    // std::cout<<"Make pairs Gen"<<std::endl;
    const unsigned int combinations = l1.size()*l2.size();
    Vec_b res(combinations);
    for (unsigned int i = 0; i < l1.size(); ++i) {
        for (unsigned int j = 0; j < l2.size(); ++j) {
            res[i*l2.size()+j] = l1[i] & l2[j];
        }
    }
    // std::cout<<"End Make pairs Gen"<<std::endl;
    return res;
}

Vec_b makePairsGen(const Vec_b& l1, const bool l2) {
    // std::cout<<"Make pairs Gen"<<std::endl;
    Vec_b res(l1.size());
    for (unsigned int i = 0; i < l1.size(); ++i) {
        res[i] = l1[i] & l2;
    }
    // std::cout<<"End Make pairs Gen"<<std::endl;
    return res;
}

bool makePairsGen(const bool c1, const bool c2) {
    return c1 & c2;
}

// --- pair muons and compute the sum of the variable
Vec_f makePairsSum(const Vec_f& l1) {
    if(l1.size() <= 1) {
        Vec_f res(0);
        return res;
    }
    const unsigned int combinations = l1.size()*(l1.size()-1)/2;
    Vec_f res(combinations);
    for (unsigned int i = 0, k=0; i < l1.size(); ++i) {
        for (unsigned int j = i+1; j < l1.size(); ++j) {
            res[k] = l1[i] + l1[j];
            k++;
        }
    }
    return res;
}

Vec_f makePairsSum(const Vec_f& l1, const Vec_f& l2) {
    // std::cout<<"End Make pairs sum"<<std::endl;
    const unsigned int combinations = l1.size()*l2.size();
    Vec_f res(combinations);
    for (unsigned int i = 0; i < l1.size(); ++i) {
        for (unsigned int j = 0; j < l2.size(); ++j) {
            res[i*l2.size()+j] = l1[i] + l2[j];
        }
    }
    // std::cout<<"End Make pairs sum"<<std::endl;
    return res;
}

Vec_f makePairsSum(const Vec_f& l1, const float l2) {
    Vec_f res(l1.size());
    for (unsigned int i = 0; i < l1.size(); ++i) {
        res[i] = l1[i] + l2;
    }
    return res;
}

float makePairsSum(const float c1, const float c2) {
    return c1 + c2;
}

}


#endif