BR_TAUToMU = 0.1739
BR_TAUToE = 0.1782
xsec_Zmm = 2001.9
xsec_Wpmunu = 11765.9
xsec_Wmmunu = 8703.87
Z_TAU_TO_LEP_RATIO = (1.-(1. - BR_TAUToMU - BR_TAUToE)**2)

xsec_13TeV = {
    'Zmumu' : xsec_Zmm,
    'Ztautau' : xsec_Zmm*Z_TAU_TO_LEP_RATIO,
    'Wplusmunu' : xsec_Wpmunu,
    'Wminusmunu' : xsec_Wmmunu,
    'Wplustaunu' : BR_TAUToMU*xsec_Wpmunu,    
    'Wminustaunu' : BR_TAUToMU*xsec_Wmmunu,
    'TTLeptonic' : 88.29,
    'TTSemileptonic' : 365.34,
    'SingleTschanLepDecays' : 3.74,
    'SingleTtWAntitop' : 19.55,
    'SingleTtchanAntitop' : 70.79,
    'SingleTtchanTop' : 119.71,
    'WW' : 75.8,
    'WZ' : 27.6,
    'ZZ2l2nu' : 0.564,
    'QCDmuEnrichPt15' : 238800,
}
