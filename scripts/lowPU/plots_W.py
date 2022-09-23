
import sys,array,math,os,copy,fnmatch
from collections import OrderedDict

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import plotter
import functions

import lz4.frame
import pickle
import narf
import hist

import hist
import numpy as np

from wremnants.datasets.datagroupsLowPU import datagroupsLowPU_Z



def doOverlow(h):

    n = h.GetNbinsX()
    h.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
    h.SetBinContent(n, h.GetBinContent(n+1) + h.GetBinContent(n))
    h.SetBinError(1, math.hypot(h.GetBinError(0), h.GetBinError(1)))
    h.SetBinError(n, math.hypot(h.GetBinError(n+1), h.GetBinError(n)))
    h.SetBinContent(0, 0)
    h.SetBinContent(n+1, 0)
    h.SetBinContent(0, 0)
    h.SetBinContent(n+1, 0)
    
    return h
    
def getSyst(In, systName, yBins):

    pass
    
def parseHist(h, norm, rebin, projectionx=[]):

    if h.InheritsFrom("TH2") and projectionx != []: h = h.ProjectionX("", projectionx[0], projectionx[1])

    if isinstance(rebin, int): h.Rebin(rebin)
    else: h = Rebin(h, rebin)
    
    h.Scale(norm)
    h = doOverlow(h)
    
    return h
    
    
def computeEGammaScaleSmearingUnc(hName, proc, fIn, fNameIn):

    systName = "egamma_scale_smear"
    
    histName_MC = "%s_%s" % (hName, proc)
    hNom_MC = copy.deepcopy(fIn.Get(histName_MC))
 
    histName_DATA = "%s_%s" % (hName, "singleMuon")
    hNom_DATA = copy.deepcopy(fIn.Get(histName_DATA))
 
    hSigmaUp = copy.deepcopy(fIn.Get("%s_syst_%sUp" % (histName_MC, systName)))
    hSigmaDown = copy.deepcopy(fIn.Get("%s_syst_%sDown" % (histName_MC, systName)))
    
    hScaleUp_DATA = copy.deepcopy(fIn.Get("%s_syst_%sUp" % (histName_DATA, systName)))
    hScaleDown_DATA = copy.deepcopy(fIn.Get("%s_syst_%sDown" % (histName_DATA, systName)))
    
    
    # project DATA differences on MC
    hScaleUp = hNom_MC.Clone("hScaleUp")
    hScaleDown = hNom_MC.Clone("hScaleDown")
    for i in range(0, hScaleUp.GetNbinsX()+2):
    
        dN = hNom_DATA.GetBinContent(i) - hScaleUp_DATA.GetBinContent(i)
        hScaleUp.SetBinContent(i, hScaleUp.GetBinContent(i) + dN)
        
        dN = hNom_DATA.GetBinContent(i) - hScaleDown_DATA.GetBinContent(i)
        hScaleDown.SetBinContent(i, hScaleDown.GetBinContent(i) + dN)
    
    
    
    fIn.cd()
    
    
    
    ## do some plotting of the ratio
    
    if hName == "m_mumu":
    
        cfg = {

            'logy'              : False,
            'logx'              : False,
            
            'xmin'              : 60,
            'xmax'              : 120,
            'ymin'              : 1e1,
            'ymax'              : 60000, # 3e6
                
            'xtitle'            : "m(#mu^{#plus},#mu^{#minus}) (GeV)",
            'ytitle'            : "Events",
                
            'topRight'          : "199 pb^{#minus1} (13 TeV)", 
            'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

            'ratiofraction'     : 0.3,
            'ytitleR'           : "Ratio",
            
            'yminR'             : 0.95,
            'ymaxR'             : 1.05,
        }
        
        
        leg = ROOT.TLegend(.20, 0.85-(len(bkgs)+2)*0.055, .5, .85)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.040)
        
        leg.AddEntry(hNom_MC, "Nominal", "L")
        leg.AddEntry(hScaleUp, "Scale up", "L")
        leg.AddEntry(hScaleDown, "Scale down", "L")

        plotter.cfg = cfg
        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio()
            
        ## top panel
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        dummyT.Draw("HIST")
           
        hNom_MC.SetLineColor(ROOT.kBlack)
        hNom_MC.SetLineWidth(2)
        hNom_MC.Draw("HIST SAME")
        
        hScaleUp.SetLineColor(ROOT.kRed)
        hScaleUp.SetLineWidth(2)
        hScaleUp.Draw("HIST SAME")
        
        hScaleDown.SetLineColor(ROOT.kBlue)
        hScaleDown.SetLineWidth(2)
        hScaleDown.Draw("HIST SAME")
        
        leg.Draw("SAME")
        
        plotter.auxRatio()  
        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        ROOT.gPad.RedrawAxis()  
        


        ## bottom panel
        canvas.cd()
        padB.Draw()
        padB.SetFillStyle(0)
        padB.cd()
        dummyB.Draw("HIST")
        dummyL.Draw("SAME")
        
        hScaleUp_ratio = hScaleUp.Clone("tmp1")
        hScaleDown_ratio = hScaleDown.Clone("tmp1")
        
        hScaleUp_ratio.Divide(hNom_MC)
        hScaleDown_ratio.Divide(hNom_MC)
         
        hScaleUp_ratio.SetLineColor(ROOT.kRed)
        hScaleUp_ratio.SetLineWidth(2)
        
        hScaleDown_ratio.SetLineColor(ROOT.kBlue)
        hScaleDown_ratio.SetLineWidth(2)

        hScaleUp_ratio.Draw("HIST SAME")
        hScaleDown_ratio.Draw("HIST SAME")

        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        ROOT.gPad.RedrawAxis()

        canvas.SaveAs("/eos/user/j/jaeyserm/www/wmass/lowPU/Zee_systs/scale.png")
        canvas.SaveAs("/eos/user/j/jaeyserm/www/wmass/lowPU/Zee_systs/scale.pdf")
        canvas.Close()    
        
    
    
    
        ##############################
        leg = ROOT.TLegend(.20, 0.85-(len(bkgs)+2)*0.055, .5, .85)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.040)
        
        leg.AddEntry(hNom_MC, "Nominal", "L")
        leg.AddEntry(hSigmaUp, "Smearing up", "L")
        leg.AddEntry(hSigmaDown, "Smearing down", "L")

        plotter.cfg = cfg
        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio()
            
        ## top panel
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        dummyT.Draw("HIST")
           
        hNom_MC.SetLineColor(ROOT.kBlack)
        hNom_MC.SetLineWidth(2)
        hNom_MC.Draw("HIST SAME")
        
        hSigmaUp.SetLineColor(ROOT.kRed)
        hSigmaUp.SetLineWidth(2)
        hSigmaUp.Draw("HIST SAME")
        
        hSigmaDown.SetLineColor(ROOT.kBlue)
        hSigmaDown.SetLineWidth(2)
        hSigmaDown.Draw("HIST SAME")
        
        leg.Draw("SAME")
        
        plotter.auxRatio()  
        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        ROOT.gPad.RedrawAxis()  
        


        ## bottom panel
        canvas.cd()
        padB.Draw()
        padB.SetFillStyle(0)
        padB.cd()
        dummyB.Draw("HIST")
        dummyL.Draw("SAME")
        
        hSigmaUp_ratio = hSigmaUp.Clone("tmp1")
        hSigmaDown_ratio = hSigmaDown.Clone("tmp1")
        
        hSigmaUp_ratio.Divide(hNom_MC)
        hSigmaDown_ratio.Divide(hNom_MC)
         
        hSigmaUp_ratio.SetLineColor(ROOT.kRed)
        hSigmaUp_ratio.SetLineWidth(2)
        
        hSigmaDown_ratio.SetLineColor(ROOT.kBlue)
        hSigmaDown_ratio.SetLineWidth(2)

        hSigmaUp_ratio.Draw("HIST SAME")
        hSigmaDown_ratio.Draw("HIST SAME")

        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        ROOT.gPad.RedrawAxis()

        canvas.SaveAs("/eos/user/j/jaeyserm/www/wmass/lowPU/Zee_systs/smearing.png")
        canvas.SaveAs("/eos/user/j/jaeyserm/www/wmass/lowPU/Zee_systs/smearing.pdf")
        canvas.Close()   
  
    
    
    
    return hScaleUp, hScaleDown, hSigmaUp, hSigmaDown


def addUnc(hTarget, hNom, hUp, hDw=None):

    for k in range(1, hTarget.GetNbinsX()+1):

        if hDw != None: sigma = 0.5*abs(abs(hUp.GetBinContent(k)-hNom.GetBinContent(k)) + abs(hDw.GetBinContent(k)-hNom.GetBinContent(k)))
        else: sigma = abs(hUp.GetBinContent(k)-hNom.GetBinContent(k))
        hTarget.SetBinError(k, math.sqrt(sigma*sigma + hTarget.GetBinError(k)*hTarget.GetBinError(k)))
                            


def recoilStatUncertainty(systName, sysHists, rebin):

    systHistName = systs[systName]['histName']
                    
    # get the nominal (only for DY)
    thn_base = fIn.Get("recoil_corr_magn_DYmumu_MiNNLO")
    h_base = thn_base.Projection(0, "E")
                    
    thn_nom = fIn.Get("RecoilSystStat_nom_DYmumu_MiNNLO") # 
    thn_pert = fIn.Get("%s_DYmumu_MiNNLO" % systHistName)
    
    groupIdx = groupSysts.index("total")        

    for iRecoilVar in range(1, 74):
                    
        thn_nom.GetAxis(1).SetRange(iRecoilVar, iRecoilVar)
        h_nom = thn_nom.Projection(0, "E")
        h_nom.SetTitle("h_nom%d" % iRecoilVar)
        h_nom.SetBinContent(0, 0) # remove underflow

        thn_pert.GetAxis(1).SetRange(iRecoilVar, iRecoilVar)
        h_pert = thn_pert.Projection(0, "E")
        h_pert.SetTitle("h_pert%d" % iRecoilVar)
        h_pert.SetBinContent(0, 0) # remove underflow
                        
        print(h_base.GetNbinsX(), h_nom.GetNbinsX(), h_pert.GetNbinsX())
                        
        hUp = copy.deepcopy(h_base)
        hUp.SetName("RecoilSystStat_para_data_m1_qTbin%dUp" % iRecoilVar)
        hUp.Add(h_nom, -1)
        hUp.Add(h_pert)
        
        hDw = copy.deepcopy(hUp)
        hDw.SetName("RecoilSystStat_para_data_m1_qTbin%dDown" % iRecoilVar)
                        
        for k in range(1, hUp.GetNbinsX()+1):
                        
            hDw.SetBinContent(k, 2*h_base.GetBinContent(k) - hUp.GetBinContent(k))
            #print(hUp.GetBinContent(k), h_base.GetBinContent(k), hUp.GetBinContent(k)/h_base.GetBinContent(k))
                        
                        
                        
        print(h_base.Integral(), h_nom.Integral(), h_pert.Integral(), h_nom.Integral()/h_pert.Integral(), hUp.Integral())
                    
        hUp = parseHist(hUp, lumi, rebin)
        hDw = parseHist(hDw, lumi, rebin)

        addUnc(sysHists[groupIdx], h_base, hUp, hDw)  
        

def parseProc(histCfg, procName, syst="", rebin=1):

    axis = histCfg['axis']
    hName = histCfg['name']
    charge = histCfg['charge']
    
    label = "%s_%s" % (hName, procName)
    groups.setHists(hName, "", label=label, procsToRead=[procName], selectSignal=False)
    bhist = groups.groups[procName][label]
    
    s = hist.tag.Slicer()
    if charge == "combined": bhist = bhist[{"charge" : s[::hist.sum]}]
    elif charge == "plus": bhist = bhist[{"charge" : bhist.axes["charge"].index(+1)}]
    elif charge == "minus": bhist = bhist[{"charge" : bhist.axes["charge"].index(-1)}]
    bhist = bhist[{"eta" : s[::hist.sum]}]
    print(bhist)
    
    rhist = narf.hist_to_root(bhist)
    rhist = functions.Rebin(rhist, rebin)
    rhist.SetName(label)
    rhist = doOverlow(rhist)
    
    if procName == "SingleMuon" or procName == "SingleElectron": pass
    else: rhist.Scale(MC_SF)
    
    #rhist.Scale(1, "width")

    print("Get histogram %s, yield=%.2f" % (label, rhist.Integral()))
    return rhist
 
 
def doRecoilStatSystParametric(histCfg, procName, hNom, h_syst, rebin):

    #return
    axis = histCfg['axis']
    hName = histCfg['name']
    import hist
    import numpy as np

    if procName not in ["WJetsToMuNu"]: return hNom, hNom
    

    tags = ["target_para", "target_perp", "source_para", "source_perp"]
    #tags = ["source_para"]
    #tags = ["target_para", "target_perp"]
    if hName == "recoil_corr_rec_para_qT": basname_syst = "recoil_para_qT_recoilSyst"
    elif hName == "recoil_corr_rec_para": basname_syst = "recoil_para_recoilSyst"
    elif hName == "recoil_corr_rec_perp": basname_syst = "recoil_perp_recoilSyst"
    elif hName == "MET_corr_rec_pt": basname_syst = "MET_recoilSyst"
    elif hName == "mT_corr_rec": basname_syst = "mT_corr_rec_recoilSyst"
    elif hName == "recoil_corr_rec_magn": basname_syst = "recoil_magn_recoilSyst"
    else: return
    
    #sys.exit("lol")
    
    s = hist.tag.Slicer()
        
    
    groups.setHists(hName, "", label=hName, procsToRead=[procName], selectSignal=False)
    bhist_nom = groups.groups[procName][hName]
    bhist_nom = bhist_nom[{"charge" : s[::hist.sum], "eta" : s[::hist.sum]}]
    rhist_nom = narf.hist_to_root(bhist_nom)
    for tag in tags:
    
        label = "%s_%s" % (basname_syst, tag)
        groups.setHists(label, "", label=label, procsToRead=[procName], selectSignal=False)
        bhist_pert = groups.groups[procName][label]
        bhist_pert = bhist_pert[{"charge" : s[::hist.sum], "eta" : s[::hist.sum]}]
        
        if procName == "SingleMuon" or procName == "SingleElectron": pass
        else: bhist_pert *= MC_SF
        
        rhist_pert = narf.hist_to_root(bhist_pert)
        for i in range(1, rhist_pert.GetNbinsY()+1):
            hUp = rhist_pert.ProjectionX("upvar%d" % i, i, i)
            hUp = functions.Rebin(hUp, rebin)
            hUp = doOverlow(hUp)
            addUnc(h_syst, hNom, hUp)
           
     
    return     
 

def doRecoilStatSyst(histCfg, procName, hNom, h_syst, rebin):

    axis = histCfg['axis']
    hName = histCfg['name']
    if procName not in ["WJetsToMuNu"]: return hNom, hNom
    
    if recoilParametric:
    
        tags = ["target_para", "target_perp", "source_para", "source_perp"]
        tags = ["source_para"]
        if hName == "recoil_corr_rec_para_qT": basname_syst = "recoil_para_qT_recoilSyst"
        elif hName == "recoil_corr_rec_para": basname_syst = "recoil_para_recoilSyst"
        elif hName == "recoil_corr_rec_perp": basname_syst = "recoil_perp_recoilSyst"
        elif hName == "MET_corr_rec_pt": basname_syst = "MET_recoilSyst"
        elif hName == "mT_corr_rec": basname_syst = "mT_corr_rec_recoilSyst"
        elif hName == "recoil_corr_rec_magn": basname_syst = "recoil_magn_recoilSyst"
        else: return
        
        s = hist.tag.Slicer()
        
        groups.setHists(hName, "", label=hName, procsToRead=[procName], selectSignal=False)
        bhist_nom = groups.groups[procName][hName]
        bhist_nom = bhist_nom[{"charge" : s[::hist.sum], "eta" : s[::hist.sum]}]
        rhist_nom = narf.hist_to_root(bhist_nom)
        for tag in tags:
            
            label = "%s_%s" % (basname_syst, tag)
            groups.setHists(label, "", label=label, procsToRead=[procName], selectSignal=False)
            bhist_pert = groups.groups[procName][label]
            bhist_pert = bhist_pert[{"charge" : s[::hist.sum], "eta" : s[::hist.sum]}]
            
            if procName == "SingleMuon" or procName == "SingleElectron": pass
            else: bhist_pert *= MC_SF
            
            rhist_pert = narf.hist_to_root(bhist_pert)
            for i in range(1, rhist_pert.GetNbinsY()+1):
                #if i != 10: continue
                hUp = rhist_pert.ProjectionX("upvar%d" % i, i, i)
                hUp = functions.Rebin(hUp, rebin)
                hUp = doOverlow(hUp)
                for j in range(0, hUp.GetNbinsX()+1):
                    print(hNom.GetBinCenter(j), hNom.GetBinContent(j), hUp.GetBinContent(j), hUp.GetBinContent(j)/hNom.GetBinContent(j) if hNom.GetBinContent(j) > 0 else 0)
                addUnc(h_syst, hNom, hUp)
                #break
    
    else:
    
        if hName == "recoil_corr_para_qT": basname_syst, axis, axis_pert = "recoil_para_qT_recoilStatUnc", "recoil_para_qT", "recoil_para_qT_pert"
        elif hName == "recoil_corr_perp": basname_syst, axis, axis_pert = "recoil_perp_recoilStatUnc", "recoil_perp", "recoil_perp_pert"
        elif hName == "MET_corr_rec_pt": basname_syst, axis, axis_pert = "MET_recoilStatUnc", "recoil_MET_pt", "recoil_MET_pt_pert"
        elif hName == "mT_corr_rec": basname_syst, axis, axis_pert = "mT_corr_rec_recoilStatUnc", "mt", "mt_pert" # recoil_corr_mt_stat_unc_
        elif hName == "recoil_corr_magn": basname_syst, axis, axis_pert = "recoil_magn_recoilStatUnc", "recoil_magn", "recoil_magn_pert"  
        else: return
         
        recoil_vars = [(1,2),(1,3),(1,4),(1,5),  (2,2),(2,3),(2,4),  (3,2),(3,3),(3,4),(3,5),  (4,2),(4,3),(4,4)]
        #recoil_vars = [(1,2)]
        #recoil_vars = []
        for k in recoil_vars:
        
            print(k)
            label = "%s_%d_%d" % (basname_syst, k[0], k[1])
            hName_syst = "%s_%d_%d" % (basname_syst, k[0], k[1])
            groups.setHists(hName_syst, "", label=label, procsToRead=[procName], selectSignal=False)
            scale_hist = groups.groups[procName][label]
            
            if procName == "SingleMuon" or procName == "SingleElectron": pass
            else: scale_hist *= MC_SF

            
            # axes: axis_recoil_para, axis_recoil_para_pert, axis_recoil_stat_unc_var
            # columns: recoil_corr_para recoil_corr_para_stat_unc_X_Y recoil_corr_stat_idx
            s = hist.tag.Slicer()
            base = scale_hist[{"recoil_stat_unc_var" : s[::hist.sum], axis_pert : s[::hist.sum], "charge" : scale_hist.axes["charge"].index(+1), "eta" : s[::hist.sum] }] # sum over qT, remove the qT axis == (RECO, M)
            nom = scale_hist[{axis_pert : s[::hist.sum], "charge" : scale_hist.axes["charge"].index(+1), "eta" : s[::hist.sum]  }] # sum over qT (perturbed)  == (RECO, M, IDX)
            pert = scale_hist[{axis : s[::hist.sum], "charge" : scale_hist.axes["charge"].index(+1), "eta" : s[::hist.sum]}] # sum over qT == (RECO_PERT IDX)
            
            #s[::hist.sum]
            #if charge == "combined": bhist = bhist[{"charge" : s[::hist.sum]}]
            #elif charge == "plus": bhist = bhist[{"charge" : bhist.axes["charge"].index(+1)}]
            #elif charge == "minus": bhist = bhist[{"charge" : bhist.axes["charge"].index(-1)}]
            
            r_base = narf.hist_to_root(base)
            r_nom = narf.hist_to_root(nom)
            r_pert = narf.hist_to_root(pert)
            
            
            print(base)
            print(nom)
            print(pert)
            
            # up histogram for each qT bin (2D histogram: param vs qT bins)
            scale_variation_hist_up = hist.Hist(*pert.axes, storage = scale_hist._storage_type(), name = hName_syst+"Up", data = pert.view(flow=True) - nom.view(flow=True) + base.view(flow=True)[..., np.newaxis])        
            rhist = narf.hist_to_root(scale_variation_hist_up)
            for i in range(1, rhist.GetNbinsY()+1):
            
                hUp = rhist.ProjectionX("upvar%d" % i, i, i)
                hUp = functions.Rebin(hUp, rebin)
                hUp = doOverlow(hUp)
                addUnc(h_syst, hNom, hUp)
        
      
 
                        
def parseHists(histCfg, leg, rebin=1, projectionx=[], noData=False):

    h_data = parseProc(histCfg, data, rebin=rebin)
    leg.AddEntry(h_data, groups.groups[data]['label'], "PE")
    
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg = h_data.Clone("bkg_nominal")
    h_bkg.Reset("ACE")
    for i,proc in enumerate(procs):
    
        hist = parseProc(histCfg, proc, rebin=rebin)
        hist.SetFillColor(ROOT.TColor.GetColor(groups.groups[proc]['color']))
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetLineStyle(1)
        
        leg.AddEntry(hist, groups.groups[proc]['label'], "F")
        st.Add(hist)
        h_bkg.Add(hist)
 
    h_bkg.SetLineColor(ROOT.kBlack)
    h_bkg.SetFillColor(0)
    h_bkg.SetLineWidth(2)

    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetLineColor(ROOT.kBlack)

    h_err = h_bkg.Clone("syst")
    h_err.SetFillColor(ROOT.kBlack)
    h_err.SetMarkerSize(0)
    h_err.SetLineWidth(0)
    h_err.SetFillStyle(3004)    
    leg.AddEntry(h_err, "Stat. + Syst. Unc.", "F")
    
    h_syst = h_bkg.Clone("syst_group%d" % i)
    h_syst.SetFillColor(ROOT.kBlack)
    h_syst.SetFillColorAlpha(ROOT.kBlack, 1)
    h_syst.SetMarkerSize(0)
    h_syst.SetLineWidth(0)
    h_syst.SetFillStyle(3004)
    for iBin in range(0, h_syst.GetNbinsX()+1): h_syst.SetBinError(iBin, 0) # reset errors for non-stat ones
            
    
    
    systs = {}
    
    # do systematics
    for i,proc in enumerate(procs):

        hNom = parseProc(histCfg, proc, rebin=rebin)
        hPert = doRecoilStatSyst(histCfg, proc, hNom, h_err, rebin)
        
        continue
        
      
        hist = None
        for proc in plotcfg[bkg]['samples']:
        
            # get nominal
            thn = None
            for hName in hNames:
                
                histName = "%s_%s" % (hName, proc)
                thn = fIn.Get(histName)
                if thn == None: continue
                else: break
        
       
            hNom = thn.Projection(0, "E")
            hNom.SetTitle("hNom_%s" % histName)
            hNom = parseHist(hNom, lumi, rebin)
            
            
            
            sys.exit()
            continue

            for systName in systs:
            
                if "RecoilSystStat" in systName or "METSystStat" in systName:

                    
                    if "METSystStat" in systName and hName == "recoil_corr_magn": continue
                    if "RecoilSystStat" in systName and hName == "MET_corr_rec_pt": continue
                    if proc != "DYmumu_MiNNLO": continue
        
                    #continue
                    systHistName = systs[systName]['histName']
                    
                    # get the nominal (only for DY)
                    if "RecoilSystStat" in systName:
                    
                        thn_nom = fIn.Get("RecoilSystStat_nom_DYmumu_MiNNLO") # binned in qT
                        thn_base = fIn.Get("recoil_corr_magn_DYmumu_MiNNLO")
                        
                    else:
                    
                        thn_nom = fIn.Get("METSystStat_nom_DYmumu_MiNNLO") # binned in qT
                        thn_base = fIn.Get("MET_corr_rec_pt_DYmumu_MiNNLO")
                     
                    h_base = thn_base.Projection(0, "E")
                    h_base = parseHist(h_base, lumi, rebin)
                     
                    # get the perturbation
                    thn_pert = fIn.Get("%s_DYmumu_MiNNLO" % systHistName)
                    
                    

                    for iRecoilVar in range(1, 40):
                        #if iRecoilVar == 40: continue
                        thn_nom.GetAxis(1).SetRange(iRecoilVar, iRecoilVar)
                        h_nom = thn_nom.Projection(0, "E")
                        h_nom.SetTitle("h_nom%d" % iRecoilVar)
                        h_nom.SetBinContent(0, 0) # remove underflow
                        h_nom.SetBinContent(h_nom.GetNbinsX()+1, 0) # remove overflow
                        h_nom = parseHist(h_nom, lumi, rebin)
                        
                        if h_nom.Integral() <= 0: continue

                        thn_pert.GetAxis(1).SetRange(iRecoilVar, iRecoilVar)
                        h_pert = thn_pert.Projection(0, "E")
                        h_pert.SetTitle("h_pert%d" % iRecoilVar)
                        h_pert.SetBinContent(0, 0) # remove underflow
                        h_pert.SetBinContent(h_pert.GetNbinsX(), 0) # remove overflow
                        h_pert.SetBinContent(h_pert.GetNbinsX()+1, 0) # remove overflow
                        h_pert = parseHist(h_pert, lumi, rebin)
                        
                        #print(h_base.GetNbinsX(), h_nom.GetNbinsX(), h_pert.GetNbinsX())
                        
                        hUp = copy.deepcopy(h_base)
                        hUp.SetName("RecoilSystStat_para_data_m1_qTbin%dUp" % iRecoilVar)
                        hUp.Add(h_nom, -1)
                        hUp.Add(h_pert)
                        
                        hDw = copy.deepcopy(hUp)
                        hDw.SetName("RecoilSystStat_para_data_m1_qTbin%dDown" % iRecoilVar)
                        
                        for k in range(1, hUp.GetNbinsX()+1):
                        
                            hDw.SetBinContent(k, 2*h_base.GetBinContent(k) - hUp.GetBinContent(k))
                            #print(hUp.GetBinCenter(k), hUp.GetBinContent(k), h_base.GetBinContent(k), hUp.GetBinContent(k)/h_base.GetBinContent(k))
                        
                        
                        #print(h_base.Integral(), h_nom.Integral(), h_pert.Integral(), hUp.Integral(), hDw.Integral())                    
                        #hUp = parseHist(hUp, lumi, rebin)
                        #hDw = parseHist(hDw, lumi, rebin)
                        
                        print(hNom.Integral(), h_base.Integral(), h_nom.Integral(), h_pert.Integral(), hUp.Integral(), hDw.Integral())
                        print(hUp.Integral()/hNom.Integral()) #### SHOULD BE 1 (as recoil does not alter the normalization)
                        
                        for kk in range(0, hUp.GetNbinsX()+1):
                            pass
                            if hNom.GetBinContent(kk) <= 0: continue
                            #print(kk, hNom.GetBinContent(kk), h_nom.GetBinContent(kk), h_pert.GetBinContent(kk), hUp.GetBinContent(kk), hUp.GetBinContent(kk)/hNom.GetBinContent(kk))
                            
                    
                        groupIdx = groupSysts.index("total")        
                        addUnc(sysHists[groupIdx], h_base, hUp, hDw)  
                    
                        #break
                    #sys.exit()
                    continue    
                
                idx = systs[systName]['idx']
                systHistName = systs[systName]['histName']
                group = systs[systName]['group']
                
                #print(systName)
                
                
                thn = None
                for hName in hNames:

                    thn = fIn.Get(histName + "_" + systHistName)
                    if thn == None: continue
                    else: break

                if len(idx) == 2: # up/down specified
                            
                    idxUp = idx[0]
                    idxDw = idx[1]
                    thn.GetAxis(1).SetRange(idxUp, idxUp)
                    hUp = thn.Projection(0, "E")
                    hUp.SetName("%s_%sUp" % (histName, systName))
                    
                    thn.GetAxis(1).SetRange(idxDw, idxDw)    
                    hDw = thn.Projection(0, "E")
                    hDw.SetName("%s_%shDw" % (histName, systName))
                    
                    hUp = parseHist(hUp, lumi, rebin)
                    hDw = parseHist(hDw, lumi, rebin)
                            
                            
                    #print(systName, bkg, proc, hNom.Integral(), hUp.Integral(), hDw.Integral())
                            
                         
                    #for k in range(1, sysHists[groupIdx].GetNbinsX()+1):

                    #    sigma = 0.5*abs(abs(hUp.GetBinContent(k)-hNom.GetBinContent(k)) + abs(hDw.GetBinContent(k)-hNom.GetBinContent(k)))
                    #    sysHists[groupIdx].SetBinError(k, math.sqrt(sigma*sigma + sysHists[groupIdx].GetBinError(k)*sysHists[groupIdx].GetBinError(k)))
                            
                else: # only 1 specified: symmetrize
                
                    idxUp = idx[0]
                    thn.GetAxis(1).SetRange(idxUp, idxUp)    
                    hUp = thn.Projection(0, "E")
                    hUp.SetName("%s_%sUp" % (histName, systName))
                    hUp = parseHist(hUp, lumi, rebin)        
                    
                    #hUp = systHist.ProjectionX("%s_%sUp" % (histName, systName), idxUp, idxUp)
                    #hUp = parseHist(hUp, lumi, rebin, projectionx=projectionx)
                    hDw = None
                    
                    #for k in range(1, hUp.GetNbinsX()+1):
                    
                    #    print(k, hNom.GetBinContent(k), hUp.GetBinContent(k), 100.*hUp.GetBinContent(k)/hNom.GetBinContent(k))
         
                    #print(systName, bkg, proc, hNom.Integral(), hUp.Integral())
                    #addUnc(sysHists[groupIdx], hNom, hUp)                    
                    #for k in range(1, sysHists[groupIdx].GetNbinsX()+1):

                    #    sigma = abs(hUp.GetBinContent(k)-hNom.GetBinContent(k))
                    #    sysHists[groupIdx].SetBinError(k, math.sqrt(sigma*sigma + sysHists[groupIdx].GetBinError(k)*sysHists[groupIdx].GetBinError(k)))
   
                
                if group in groupSysts:
                    groupIdx = groupSysts.index(group)        
                    addUnc(sysHists[groupIdx], hNom, hUp, hDw)   
                
                if "total" in groupSysts:
                    groupIdx = groupSysts.index("total")        
                    addUnc(sysHists[groupIdx], hNom, hUp, hDw)   
          

                    ## do EGammaScaleSmearing
                    '''
                    hScaleUp, hScaleDown, hSigmaUp, hSigmaDown = computeEGammaScaleSmearingUnc(hName, proc, fIn, fNameIn)
                    
                    hScaleUp = parseHist(hScaleUp, lumi, rebin)
                    hScaleDown = parseHist(hScaleDown, lumi, rebin)
                    hSigmaUp = parseHist(hSigmaUp, lumi, rebin)
                    hSigmaDown = parseHist(hSigmaDown, lumi, rebin)
                    
                    print("EGammaScale", bkg, proc, hNom.Integral(), hScaleUp.Integral(), hScaleDown.Integral())
                    print("EGammaSmear", bkg, proc, hNom.Integral(), hSigmaUp.Integral(), hSigmaDown.Integral())
                    #sys.exit()
                    for k in range(1, h_bkg_err.GetNbinsX()+1):
                        #continue
                        sigma = 0.5*abs(abs(hSigmaUp.GetBinContent(k)-hNom.GetBinContent(k)) + abs(hSigmaDown.GetBinContent(k)-hNom.GetBinContent(k)))
                        h_bkg_err.SetBinError(k, math.sqrt(sigma*sigma + h_bkg_err.GetBinError(k)*h_bkg_err.GetBinError(k)))
                    
                    for k in range(1, h_bkg_err.GetNbinsX()+1):
                        #continue
                        sigma = 0.5*abs(abs(hScaleUp.GetBinContent(k)-hNom.GetBinContent(k)) + abs(hScaleDown.GetBinContent(k)-hNom.GetBinContent(k)))
                        h_bkg_err.SetBinError(k, math.sqrt(sigma*sigma + h_bkg_err.GetBinError(k)*h_bkg_err.GetBinError(k)))
                    '''
       
    # ratios (bands)
    h_bkg_ratio = h_bkg.Clone("h_bkg_ratio") # nominal point, need to remove stat. error
    h_err_ratio = h_err.Clone("syst_ratio")
    h_err_ratio.Divide(h_bkg_ratio)


    # Data/MC ratio (the error bars represent the data uncertainty)
    h_ratio = h_data.Clone("h_ratio")
    for i in range(h_bkg .GetNbinsX()+1): h_bkg.SetBinError(i, 0) # set MC errors to zero
    h_ratio.Divide(h_bkg)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(0.7)
    h_ratio.SetMarkerColor(ROOT.kBlack)
    h_ratio.SetLineColor(ROOT.kBlack)
    h_ratio.SetLineWidth(1)
    
    return h_data, st, h_bkg, h_ratio, h_err, h_err_ratio



def singlePlot(histCfg, fOut, xMin, xMax, yMin, yMax, xLabel, yLabel, logY=True, rebin=1, legPos=[], projectionx=[], yRatio = 1.3):
    
    # default cfg
    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : 60,
        'xmax'              : 120,
        'ymin'              : 1e1,
        'ymax'              : 1e5, # 3e6
            
        'xtitle'            : "m(#mu,#mu) (GeV)",
        'ytitle'            : "Events",
            
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

        'ratiofraction'     : 0.3,
        'ytitleR'           : "Data/MC",
        
        'yminR'             : 1-(yRatio-1), # 0.7
        'ymaxR'             : yRatio, # 1.3
    }
    
    
    leg = ROOT.TLegend(.20, 0.88-(len(procs)+2)*0.05, .5, .88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.040)
    
    h_data, st, h_bkg, h_ratio, h_err, h_err_ratio = parseHists(histCfg, leg, rebin, projectionx=projectionx)

    cfg['logy'] = logY
    cfg['xmin'] = xMin
    cfg['xmax'] = xMax
    cfg['ymin'] = yMin
    cfg['ymax'] = yMax
    
    cfg['xtitle'] = xLabel
    cfg['ytitle'] = yLabel
    

    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio()
        
    ## top panel
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    dummyT.Draw("HIST")
        
    st.Draw("HIST SAME")
    
    h_err.Draw("E2 SAME")
    h_bkg.Draw("HIST SAME")
    #h_data.Draw("PE SAME")
    leg.Draw("SAME")
    
    plotter.auxRatio()  
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()  
    

    ## bottom panel
    canvas.cd()
    padB.Draw()
    padB.SetFillStyle(0)
    padB.cd()
    dummyB.Draw("HIST")
    dummyL.Draw("SAME")
    
    #h_ratio.Draw("P SAME") # E2 SAME
    h_err_ratio.Draw("E2 SAME")

    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("%s/%s.png" % (outDir, fOut))
    canvas.SaveAs("%s/%s.pdf" % (outDir, fOut))
    canvas.Close()



def singlePlot_noData(hName, fOut, xMin, xMax, yMin, yMax, xLabel, yLabel, logY=True, rebin=1, legPos=[], projectionx=[]):
    
    # default cfg
    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : 60,
        'xmax'              : 120,
        'ymin'              : 1e1,
        'ymax'              : 1e5, # 3e6
            
        'xtitle'            : "m(#mu,#mu) (GeV)",
        'ytitle'            : "Events",
            
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

    }
    
    
    leg = ROOT.TLegend(.20, 0.90-(len(bkgs)+1)*0.05, .55, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.040)
    
    h_data, st, h_bkg, sysHists, h_ratio, sysHists_ratio = parseHists(hName, leg, rebin, projectionx=projectionx, noData=True)
    #h_data, st, h_bkg, h_bkg_err, h_ratio, h_bkg_err_ratio = parseHists(hName, leg, rebin)

    cfg['logy'] = logY
    cfg['xmin'] = xMin
    cfg['xmax'] = xMax
    cfg['ymin'] = yMin
    cfg['ymax'] = yMax
    
    cfg['xtitle'] = xLabel
    cfg['ytitle'] = yLabel
    

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()
        
    
    canvas.SetGrid()
    dummy.Draw("HIST")
        
    st.Draw("HIST SAME")
    
    #for h in sysHists: h.Draw("E2 SAME")
    sysHists[0].Draw("E2 SAME")
    h_bkg.Draw("HIST SAME")
    leg.Draw("SAME")
    
    plotter.aux(canvas)  
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()  
    


    canvas.SaveAs("%s/%s.png" % (outDir, fOut))
    canvas.SaveAs("%s/%s.pdf" % (outDir, fOut))
    canvas.Close()


def MET2d():

    histCfg = {"name": "recoil_corr_rec_para_perp_2d", "axis": "recoil_perp" }
    h_data = parseProc(histCfg, data)
    hist_mc = None
    for i,proc in enumerate(procs):
    
        hist = parseProc(histCfg, proc)
        if hist_mc == None: hist_mc = hist
        else: hist_mc.Add(hist)
 
    h_data.Divide(hist_mc)
    
    for i in range(1, h_data.GetNbinsX()+1):
        for j in range(1, h_data.GetNbinsX()+1):
        
            val = h_data.GetBinContent(i, j)
            if val >= 1 and val < 1.05: newval = 1.025
            elif val >= 1.05 and val < 1.1: newval = 1.075
            elif val >= 1.1 and val < 1.15: newval = 1.125
            elif val >= 1.15 and val < 1.2: newval = 1.175
            elif val >= 1.2: newval = 1.3
            elif val <= 1 and val > 0.95: newval = 0.975
            elif val <= 0.95 and val > 0.9: newval = 0.925
            elif val <= 0.9 and val > 0.85: newval = 0.875
            elif val <= 0.85 and val > 0.8: newval = 0.825
            elif val <= 0.8: newval = 0.7
            h_data.SetBinContent(i, j, newval)
            #print("%.3f " % h_data.GetBinContent(i, j), end="")

    cfg = {

        'logy'              : False,
        'logx'              : False,
    
        'xmin'              : -10,
        'xmax'              : 10,
        'ymin'              : -10,
        'ymax'              : 10,
        
        'xtitle'            : "Recoil parallel",
        'ytitle'            : "Recoil perp",
        
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",
    }
    

    
    #ROOT.gStyle.SetPalette(6, [ROOT.kRed+2, ROOT.kOrange+2 , ROOT.kGreen +2 , ROOT.kGreen +2 , ROOT.kAzure +2, ROOT.kBlue  +2])
    #ROOT.gStyle.SetNumberContours(6)
    #ROOT.gStyle.SetPalette(ROOT.kVisibleSpectrum)
    

    plotter.cfg = cfg
    canvas = plotter.canvas()
    canvas.SetRightMargin(0.12)
    #canvas.SetLogz()
    dummy = plotter.dummy()  
    canvas.cd()
    dummy.Draw("HIST")
    
    ROOT.gStyle.SetPaintTextFormat("4.3f")
    h_data.SetMarkerSize(0.4) 
    
    h_data.GetZaxis().SetRangeUser(0.8, 1.2)
    #h_data.Draw("SAME COLZ TEXTE")
    h_data.Draw("SAME COLZ")      
    
    ells = []
    for r in [1, 2, 3]:
        met1 = ROOT.TEllipse(0., 0., r, r)
        met1.SetFillStyle(4000)
        met1.SetLineWidth(1)
        met1.Draw("SAME")
        ells.append(met1)

    
    plotter.aux(canvas)
    
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.RedrawAxis()
    canvas.SaveAs("%s/MET_2d.png" % (outDir))
    canvas.SaveAs("%s/MET_2d.pdf" % (outDir))
    canvas.Delete()
    
    for i in range(1, h_data.GetNbinsX()+1):
        for j in range(1, h_data.GetNbinsX()+1):
        
            print("%.3f " % h_data.GetBinContent(i, j), end="")
        print()
    
        #hx = h_data.ProjectionX("x", i, i)
        #hy = h_data.ProjectionY("y", i, i)
        
        #print(hx.Integral()/hx.GetNbinsX(), hy.Integral()/hy.GetNbinsX())
    
 
            
	
if __name__ == "__main__":

    print("Start")
    flavor = "mu"
    met = "RawPFMET" # DeepMETReso RawPFMET
    charge = "combined" # combined plus minus
    recoilParametric = True
    #flavor = "ee"
    
    MC_SF = 1.0
    if flavor == "mu": MC_SF = 1.026
    MC_SF = 1.0
    

    #outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/Z%s/plots_DeepMET/" % flavor
    if recoilParametric: outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/W%s/plots_%s_%s_parametric/" % (flavor, charge, met)
    else: outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/W%s/plots_%s_%s_binned/" % (flavor, charge, met)
    functions.prepareDir(outDir, remove=True)
    
    print("Open")
    groups = datagroupsLowPU_Z("lowPU_%s_%s.pkl.lz4" % (flavor, met), flavor=flavor)


    if flavor == "mu":
    
        procs = ['VV', 'DY', 'WJetsToTauNu',  'TTbar', 'WJetsToMuNu'] 
        procs = ['EWK',  'TTbar', 'WJetsToMuNu'] 
        data = 'SingleMuon'
        
        
    if flavor == "ee":
    
        procs = ['EWK', 'TTbar', 'DYee']
        data = 'SingleElectron'    
    print("Plotting") 
    
    
    #hist = parseProc({"name": "mT_corr_rec_before", "axis": "mt", "charge": charge }, "WJetsToMuNu")
    #for i in range(0, hist.GetNbinsX()+1):
    #    print(i, hist.GetBinCenter(i), hist.GetBinContent(i))
    
    #sys.exit()
    

    #singlePlot({"name": "mT_uncorr", "axis": "mt", "charge": charge }, "mT_uncorr", 40, 200, 1e0, 1e8, "m_{T} (GeV) (uncorrected)", "Events", rebin=1)
    #singlePlot({"name": "mT_corr_xy", "axis": "mt", "charge": charge }, "mT_corr_xy", 40, 200, 1e0, 1e8, "m_{T} (GeV) (xy corrected)", "Events", rebin=1)
    singlePlot({"name": "mT_corr_rec", "axis": "mt", "charge": charge }, "mT_corr_rec", 40, 200, 1e0, 1e6, "m_{T} (GeV) (recoil corrected)", "Events", rebin=1, yRatio=1.08)
    

    sys.exit()
    
    #bins_MET = [0, 0.5, 1, 1.5] + list(range(2, 20, 1)) + list(range(20, 50, 2)) + list(range(50, 70, 4)) + [70, 80, 90, 100] # was 2 before
    bins_MET = list(range(0, 10, 2)) + list(range(10, 70, 1)) + list(range(70, 100, 2)) + [100]
    singlePlot({"name": "MET_uncorr_pt", "axis": "MET_pt", "charge": charge  }, "MET_uncorr_pt", 0, 100, 1, 1e8, "MET p_{T} (uncorrected)", "Events", rebin=bins_MET)
    singlePlot({"name": "MET_corr_xy_pt", "axis": "MET_pt", "charge": charge  }, "MET_corr_xy_pt", 0, 100, 1, 1e8, "MET p_{T} (xy corrected)", "Events", rebin=bins_MET)
    singlePlot({"name": "MET_corr_rec_pt", "axis": "MET_pt", "charge": charge  }, "MET_corr_rec_pt", 0, 100, 1, 1e8, "MET p_{T} (recoil corrected)", "Events", rebin=bins_MET)
    
    sys.exit()
    
    
    #singlePlot({"name": "massZ", "axis": "mll" }, "massZ", 60, 120, 1e0, 1e6, "m(l, l) (GeV)", "Events")
    
    #singlePlot({"name": "counts", "axis": "MET_xy" }, "counts", 0, 5, 0, 150, "MET y", "Events", rebin=1, logY=False)
    #singlePlot({"name": "counts_weighted", "axis": "MET_xy" }, "counts_weighted", 0, 5, 0, 150, "MET y", "Events", rebin=1, logY=False)
 
    ##### MET PLOTS
    bins_MET = list(range(0, 50, 2)) + list(range(50, 70, 4)) + [70, 80, 90, 100] # was 2 before
    singlePlot({"name": "MET_uncorr_pt", "axis": "MET_pt" }, "MET_uncorr_pt", 0, 100, 1e1, 1e7, "MET p_{T} (uncorrected)", "Events", rebin=bins_MET)
    singlePlot({"name": "MET_corr_lep_pt", "axis": "MET_pt" }, "MET_corr_lep_pt", 0, 100, 1e1, 1e7, "MET p_{T} (lepton corrected)", "Events", rebin=bins_MET)
    singlePlot({"name": "MET_corr_xy_pt", "axis": "MET_pt" }, "MET_corr_xy_pt", 0, 100, 1e1, 1e7, "MET p_{T} (XY corrected)", "Events", rebin=bins_MET)
    #singlePlot({"name": "MET_corr_xy_pt_cut", "axis": "MET_pt" }, "MET_corr_xy_pt_cut", 0, 100, 1, 1e7, "MET p_{T} (XY corrected)", "Events", rebin=bins_MET)
    #singlePlot({"name": "MET_corr_rec_pt", "axis": "MET_pt" }, "MET_corr_rec_pt", 0, 100, 1, 1e7, "MET p_{T} (recoil corrected)", "Events", rebin=bins_MET)
    
    sys.exit()
    singlePlot({"name": "MET_uncorr_phi", "axis": "recoil_MET_phi" }, "MET_uncorr_phi", -4, 4, 1e3, 1e7, "MET #phi (uncorrected)", "Events", rebin=1)
    singlePlot({"name": "MET_corr_lep_phi", "axis": "recoil_MET_phi" }, "MET_corr_lep_phi", -4, 4, 1e3, 1e7, "MET #phi (lepton corrected)", "Events", rebin=1)
    singlePlot({"name": "MET_corr_xy_phi", "axis": "recoil_MET_phi" }, "MET_corr_xy_phi", -4, 4, 1e3, 1e7, "MET #phi (XY corrected)", "Events", rebin=1)
    #singlePlot({"name": "MET_corr_rec_phi", "axis": "recoil_MET_phi" }, "MET_corr_rec_phi", -4, 4, 1e1, 1e6, "MET #phi (recoil corrected)", "Events", rebin=1)
    #singlePlot({"name": "MET_corr_xy_phi_cut", "axis": "recoil_MET_phi" }, "MET_corr_xy_phi_cut", -4, 4, 1e1, 1e6, "MET #phi (XY corrected)", "Events", rebin=1)
 
    singlePlot({"name": "lep_iso", "axis": "iso" }, "lep_iso", 0, 2, 1e4, 1e9, "Lepton isolation", "Events", rebin=1)

 
    sys.exit()
    
    

    
    ## recoil plots
    bins_recoil_para_perp = [-10000, -100, -80, -70, -60, -50, -46, -42, -38, -34] + list(range(-30, 30, 2)) + [30, 34, 38, 42, 46, 50, 60, 70, 80, 100, 10000]
    bins_recoil_magn = list(range(0, 100, 4)) + [100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
    singlePlot({"name": "recoil_uncorr_para_qT", "axis": "recoil_perp" }, "recoil_uncorr_para_qT", -100, 100, 1e-1, 1e6, "U_{#parallel} + q_{T} (GeV) (uncorrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_uncorr_perp", "axis": "recoil_perp" }, "recoil_uncorr_perp", -100, 100, 1e-1, 1e6, "U_{#perp}  (GeV) (uncorrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_uncorr_para_qT_perp", "axis": "recoil_para_qT" }, "recoil_uncorr_para_qT_perp", -100, 100, 1e-1, 1e6, "U_{#parallel} + U_{#perp}  (GeV) (uncorrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_uncorr_magn", "axis": "recoil_magn" }, "recoil_uncorr_magn", 0, 200, 1e-1, 1e6, "|U| (GeV) (uncorrected)", "Events", rebin=bins_recoil_magn)
    
    singlePlot({"name": "recoil_corr_lep_para_qT", "axis": "recoil_perp" }, "recoil_corr_lep_para_qT", -100, 100, 1e-1, 1e6, "U_{#parallel} + q_{T} (GeV) (lepton corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_lep_perp", "axis": "recoil_perp" }, "recoil_corr_lep_perp", -100, 100, 1e-1, 1e6, "U_{#perp}  (GeV) (lepton corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_lep_para_qT_perp", "axis": "recoil_para_qT" }, "recoil_corr_lep_para_qT_perp", -100, 100, 1e-1, 1e6, "U_{#parallel} + U_{#perp}  (GeV) (lepton corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_lep_magn", "axis": "recoil_magn" }, "recoil_corr_lep_magn", 0, 200, 1e-1, 1e6, "|U| (GeV) (lepton corrected)", "Events", rebin=bins_recoil_magn)
    
    
    singlePlot({"name": "recoil_corr_xy_para_qT", "axis": "recoil_perp" }, "recoil_corr_xy_para_qT", -100, 100, 1e-1, 1e6, "U_{#parallel} + q_{T} (GeV) (XY corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_xy_perp", "axis": "recoil_perp" }, "recoil_corr_xy_perp", -100, 100, 1e-1, 1e6, "U_{#perp}  (GeV) (XY corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_xy_para_qT_perp", "axis": "recoil_para_qT" }, "recoil_corr_xy_para_qT_perp", -100, 100, 1e-1, 1e6, "U_{#parallel} + U_{#perp}  (GeV) (XY corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_xy_magn", "axis": "recoil_magn" }, "recoil_corr_xy_magn", 0, 200, 1e-1, 1e6, "|U| (GeV) (XY corrected)", "Events", rebin=bins_recoil_magn)
    
    singlePlot({"name": "recoil_corr_rec_para_qT", "axis": "recoil_para_qT" }, "recoil_corr_rec_para_qT", -100, 100, 1e-1, 1e6, "U_{#parallel} + q_{T} (GeV) (recoil corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_rec_perp", "axis": "recoil_perp" }, "recoil_corr_rec_perp", -100, 100, 1e-1, 1e6, "U_{#perp}  (GeV) (recoil corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_rec_para_qT_perp", "axis": "recoil_para_qT" }, "recoil_corr_rec_para_qT_perp", -100, 100, 1e-1, 1e6, "U_{#parallel} + U_{#perp}  (GeV) (recoil corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_rec_magn", "axis": "recoil_magn" }, "recoil_corr_rec_magn", 0, 200, 1e-1, 1e6, "|U| (GeV) (recoil corrected)", "Events", rebin=bins_recoil_magn)
    
 
    singlePlot({"name": "recoil_corr_xy_para_qT_cut", "axis": "recoil_perp" }, "recoil_corr_xy_para_qT_cut", -100, 100, 1e-1, 1e6, "U_{#parallel} + q_{T} (GeV) (XY corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_xy_perp_cut", "axis": "recoil_perp" }, "recoil_corr_xy_perp_cut", -100, 100, 1e-1, 1e6, "U_{#perp}  (GeV) (XY corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_xy_para_qT_perp_cut", "axis": "recoil_para_qT" }, "recoil_corr_xy_para_qT_perp_cut", -100, 100, 1e-1, 1e6, "U_{#parallel} + U_{#perp}  (GeV) (XY corrected)", "Events", rebin=bins_recoil_para_perp)
    singlePlot({"name": "recoil_corr_xy_magn_cut", "axis": "recoil_magn" }, "recoil_corr_xy_magn_cut", 0, 200, 1e-1, 1e6, "|U| (GeV) (XY corrected)", "Events", rebin=bins_recoil_magn)
    
    singlePlot({"name": "recoil_corr_xy_para_cut", "axis": "recoil_para" }, "recoil_corr_xy_para_cut", -200, 100, 1e-1, 1e6, "U_{#parallel} (GeV)", "Events", rebin=5)
 
    
    singlePlot({"name": "recoil_corr_xy_para", "axis": "recoil_para" }, "recoil_corr_xy_para", -200, 100, 1e-1, 1e6, "U_{#parallel} (GeV)", "Events", rebin=5)
    singlePlot({"name": "recoil_corr_rec_para", "axis": "recoil_para" }, "recoil_corr_rec_para", -200, 100, 1e-1, 1e6, "U_{#parallel} (GeV)", "Events", rebin=5)
    
    bins_recoil_para_perp_abs = list(range(0, 10, 1)) + list(range(10, 30, 2)) + [30, 34, 38, 42, 46, 50, 60, 70, 80, 100, 10000]
    singlePlot({"name": "recoil_corr_rec_para_qT_abs", "axis": "recoil_para_qT" }, "recoil_corr_rec_para_qT_abs", 0, 100, 1e-1, 1e6, "U_{#parallel} + q_{T} (GeV) (recoil corrected)", "Events", rebin=bins_recoil_para_perp_abs)
    singlePlot({"name": "recoil_corr_rec_perp_abs", "axis": "recoil_perp" }, "recoil_corr_rec_perp_abs", 0, 100, 1e-1, 1e6, "U_{#perp}  (GeV) (recoil corrected)", "Events", rebin=bins_recoil_para_perp_abs)

    singlePlot({"name": "recoil_corr_rec_para_plus_qt_squared", "axis": "recoil_para_qT" }, "recoil_corr_rec_para_plus_qt_squared", 0, 100, 1e-1, 1e6, "U_{#parallel} + q_{T} (GeV) (recoil corrected)", "Events")
    
    singlePlot({"name": "recoil_corr_rec_para_qT_squared", "axis": "recoil_para_qT" }, "recoil_corr_rec_para_qT_squared", 0, 100, 1e-1, 1e6, "U_{#parallel} + q_{T} (GeV) (recoil corrected)", "Events")
    singlePlot({"name": "recoil_corr_rec_perp_squared", "axis": "recoil_perp" }, "recoil_corr_rec_perp_squared", 0, 100, 1e-1, 1e6, "U_{#perp}  (GeV) (recoil corrected)", "Events")
    

    

    
    mT_bins = [0, 10, 15, 20, 25, 30, 35,] + list(range(40, 100, 2)) + [100, 102, 104, 106, 108, 110, 115, 120, 125, 130, 140, 160, 200]
    singlePlot({"name": "mT_uncorr", "axis": "mt" }, "mT_uncorr", 0, 200, 1e-1, 1e6, "m_{T} (GeV) (uncorrected)", "Events", rebin=mT_bins)
    singlePlot({"name": "mT_corr_lep", "axis": "mt" }, "mT_corr_lep", 0, 200, 1e-1, 1e6, "m_{T} (GeV) (lepton corrected)", "Events", rebin=mT_bins)
    singlePlot({"name": "mT_corr_xy", "axis": "mt" }, "mT_corr_xy", 0, 200, 1e-1, 1e6, "m_{T} (GeV) (XY corrected)", "Events", rebin=mT_bins)
    singlePlot({"name": "mT_corr_rec", "axis": "mt" }, "mT_corr_rec", 0, 200, 1e-1, 1e6, "m_{T} (GeV)  (recoil corrected)", "Events", rebin=mT_bins)
    
    #singlePlot({"name": "recoil_uncorr_para", "axis": "recoil_para" }, "recoil_uncorr_para", -200, 100, 1e-1, 1e6, "Recoil (parallel)", "Events", rebin=recoil_bins)
    
    bins = [-200, -100, -80, -70, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 200]
    
    
    bins = [-200, -100, -75, -50, -40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100, 200]
    
    recoil_qTbins = list(range(0, 50, 1)) + list(range(50, 70, 2)) + list(range(70, 100, 5)) + list(range(100, 150, 10)) + [150, 200, 300]
    singlePlot({"name": "qT", "axis": "qT" }, "qT", 0, 300, 0.1, 1e6, "q_{T} (GeV)", "Events", rebin=recoil_qTbins)
    singlePlot({"name": "qT_qTreweighted", "axis": "qT" }, "qT_qTreweighted", 0, 100, 0.1, 1e6, "q_{T} (GeV)", "Events", rebin=recoil_qTbins, yRatio=1.05)
 
    singlePlot({"name": "METx_uncorr", "axis": "MET_xy" }, "METx_uncorr", -100, 100, 1e1, 1e6, "MET x (uncorrected)", "Events", rebin=1)
    singlePlot({"name": "METy_uncorr", "axis": "MET_xy" }, "METy_uncorr", -100, 100, 1e1, 1e6, "MET y (uncorrected)", "Events", rebin=1)
    
    singlePlot({"name": "METx_corr_lep", "axis": "MET_xy" }, "METx_corr_lep", -100, 100, 1e1, 1e6, "MET x (lepton corrected)", "Events", rebin=1)
    singlePlot({"name": "METy_corr_lep", "axis": "MET_xy" }, "METy_corr_lep", -100, 100, 1e1, 1e6, "MET y (lepton corrected)", "Events", rebin=1)
 
    singlePlot({"name": "METx_corr_xy", "axis": "MET_xy" }, "METx_corr_xy", -100, 100, 1e1, 1e6, "MET x (XY corrected)", "Events", rebin=1)
    singlePlot({"name": "METy_corr_xy", "axis": "MET_xy" }, "METy_corr_xy", -100, 100, 1e1, 1e6, "MET y (XY corrected)", "Events", rebin=1)
    
    singlePlot({"name": "METx_corr_rec", "axis": "MET_xy" }, "METx_corr_rec", -100, 100, 1e1, 1e6, "MET x (recoil corrected)", "Events", rebin=1)
    singlePlot({"name": "METy_corr_rec", "axis": "MET_xy" }, "METy_corr_rec", -100, 100, 1e1, 1e6, "MET y (recoil corrected)", "Events", rebin=1)
    
    singlePlot({"name": "npv", "axis": "recoil_npv" }, "npv", 0, 15, 1e1, 1e7, "NPV", "Events", rebin=1, yRatio=1.8)
 
 
    singlePlot({"name": "MET_corr_xy_ll_dPhi", "axis": "recoil_MET_phi" }, "MET_corr_xy_ll_dPhi", -4, 4, 1e1, 1e6, "#Delta#phi(MET xy corr, ll) (lepton corrected)", "Events", rebin=1)
    singlePlot({"name": "MET_corr_rec_ll_dPhi", "axis": "recoil_MET_phi" }, "MET_corr_rec_ll_dPhi", -4, 4, 1e1, 1e6, "#Delta#phi(MET rec corr, ll) (lepton corrected)", "Events", rebin=1)
    singlePlot({"name": "ll_phi", "axis": "recoil_MET_phi" }, "ll_phi", -4, 4, 1e1, 1e6, "ll #phi (lepton corrected)", "Events", rebin=1)
    singlePlot({"name": "MET_corr_rec_xy_dPhi", "axis": "recoil_MET_phi" }, "MET_corr_rec_xy_dPhi", -4, 4, 1e1, 1e6, "#phi(MET rec corr) - #phi(MET xy corr)", "Events", rebin=1)
 
    singlePlot({"name": "counts", "axis": "MET_xy" }, "counts", 0, 5, 0, 150, "MET y", "Events", rebin=1, logY=False)
    singlePlot({"name": "counts_weighted", "axis": "MET_xy" }, "counts_weighted", 0, 5, 0, 150, "MET y", "Events", rebin=1, logY=False)
 
    sys.exit()
    
    


    singlePlot({"name": "METx_rec_xy_corr", "axis": "MET_xy" }, "METx_rec_xy_corr", -100, 100, 1e1, 1e6, "MET x", "Events", rebin=1)
    singlePlot({"name": "METy_rec_xy_corr", "axis": "MET_xy" }, "METy_rec_xy_corr", -100, 100, 1e1, 1e6, "MET y", "Events", rebin=1)

    singlePlot({"name": "MET_corr_rec_xy_pt", "axis": "MET_pt" }, "MET_corr_rec_xy_pt", 0, 100, 1, 1e7, "#phi MET", "Events", rebin=bins_MET)
    singlePlot({"name": "MET_corr_rec_xy_phi", "axis": "recoil_MET_phi" }, "MET_corr_rec_xy_phi", -4, 4, 1e1, 1e6, "#phi MET", "Events", rebin=1)
    

    
    


    
    sys.exit()
    
    singlePlot({"name": "massZ", "axis": "mll" }, "massZ_nolog", 60, 120, 0, 2.5e4, "m(l, l) (GeV)", "Events", logY=False)
    singlePlot({"name": "massZ", "axis": "mll" }, "massZ", 60, 120, 1e0, 1e6, "m(l, l) (GeV)", "Events")
    mZ_bins = [60, 65, 70, 72, 74, 76, 78, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 102, 104, 106, 108, 110, 115, 120]
    singlePlot({"name": "massZ", "axis": "mll" }, "massZ_binned", 60, 120, 1e0, 1e6, "m(l, l) (GeV)", "Events", rebin=mZ_bins, yRatio = 1.15)
    
    
    
    singlePlot({"name": "yZ", "axis": "yll" }, "yZ", -2.5, 2.5, 1e1, 1e6, "Rapidity (l, l)", "Events", rebin=1)
    singlePlot({"name": "ptZ", "axis": "ptll" }, "ptZ", 0, 300, 1e0, 1e5, "p_{T}(l, l) (GeV)", "Events", rebin=1)
    
    singlePlot({"name": "lep_pt", "axis": "ptl" }, "lep_pt", 0, 200, 1e0, 1e7, "p_{T} leptons (GeV)", "Events", rebin=1)
    singlePlot({"name": "lep_pt_leading", "axis": "ptl" }, "lep_pt_leading", 0, 200, 1e0, 1e7, "p_{T} leading lepton (GeV)", "Events", rebin=1)
    singlePlot({"name": "lep_pt_subleading", "axis": "ptl" }, "lep_pt_subleading", 0, 200, 1e0, 1e7, "p_{T} subleading lepton (GeV)", "Events", rebin=1)
    singlePlot({"name": "lep_pt_plus", "axis": "ptl" }, "lep_pt_plus", 0, 200, 1e0, 1e7, "p_{T} positive lepton (GeV)", "Events", rebin=1)
    singlePlot({"name": "lep_pt_minus", "axis": "ptl" }, "lep_pt_minus", 0, 200, 1e0, 1e7, "p_{T} negative lepton (GeV)", "Events", rebin=1)
    singlePlot({"name": "lep_pt_trg", "axis": "ptl" }, "lep_pt_trg", 0, 200, 1e0, 1e7, "p_{T} triggered lepton(s)  (GeV)", "Events", rebin=1)
    singlePlot({"name": "lep_pt_nontrg", "axis": "ptl" }, "lep_pt_nontrg", 0, 200, 1e0, 1e7, "p_{T} non-triggered lepton (GeV)", "Events", rebin=1)
    
    singlePlot({"name": "lep_eta", "axis": "etal" }, "lep_eta", -2.5, 2.5, 1e1, 1e6, "#eta leptons", "Events", rebin=1)
    singlePlot({"name": "lep_eta_leading", "axis": "etal" }, "lep_eta_leading", -2.5, 2.5, 1e1, 1e6, "#eta leading lepton", "Events", rebin=1)
    singlePlot({"name": "lep_eta_subleading", "axis": "etal" }, "lep_eta_subleading", -2.5, 2.5, 1e1, 1e6, "#eta subleading lepton", "Events", rebin=1)
    singlePlot({"name": "lep_eta_plus", "axis": "etal" }, "lep_eta_plus", -2.5, 2.5, 1e1, 1e6, "#eta positive lepton", "Events", rebin=1)
    singlePlot({"name": "lep_eta_minus", "axis": "etal" }, "lep_eta_minus", -2.5, 2.5, 1e1, 1e6, "#eta negative lepton", "Events", rebin=1)
    singlePlot({"name": "lep_eta_trg", "axis": "etal" }, "lep_eta_trg", -2.5, 2.5, 1e1, 1e6, "#eta triggered lepton(s)", "Events", rebin=1)
    singlePlot({"name": "lep_eta_nontrg", "axis": "etal" }, "lep_eta_nontrg", -2.5, 2.5, 1e1, 1e6, "#eta non-triggered lepton", "Events", rebin=1)
    
    singlePlot({"name": "noTrigMatch", "axis": "lin" }, "noTrigMatch", 0, 5, 100, 1e7, "Number of matched triggered objects", "Events", rebin=1)
    
    
    singlePlot({"name": "qT", "axis": "recoil_qT" }, "qT", 0, 300, 0.1, 1e6, "q_{T} (GeV)", "Events", rebin=1)
    
    

    
    singlePlot({"name": "RawMET_sumEt", "axis": "RawMET_sumEt" }, "RawMET_sumEt", 0, 1000, 1, 1e7, "#Sigma E_{T}", "Events", rebin=5)
    
    
    singlePlot({"name": "recoil_corr_magn", "axis": "recoil_magn" }, "recoil_corr_magn", 0, 300, 1e-1, 1e6, "|U| (GeV)", "Events", rebin=5)
    singlePlot({"name": "recoil_corr_para", "axis": "recoil_para" }, "recoil_corr_para", -200, 100, 1e-1, 1e6, "U_{#parallel} (GeV)", "Events", rebin=5)
    bins = [-200, -100, -80, -70, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 200]
    singlePlot({"name": "recoil_corr_para_qT", "axis": "recoil_para_qT" }, "recoil_corr_para_qT", -100, 100, 1e-1, 1e6, "U_{#parallel} + q_{T} (GeV)", "Events", rebin=bins)
    #singlePlot({"name": "recoil_corr_perp", "axis": "recoil_perp" }, "recoil_corr_perp", -200, 200, 1e-1, 1e6, "U_{#perp}  (GeV)", "Events", rebin=5)
    
    bins = [-200, -100, -75, -50, -40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100, 200]
    singlePlot({"name": "recoil_corr_perp", "axis": "recoil_perp" }, "recoil_corr_perp", -100, 100, 1e-1, 1e6, "U_{#perp}  (GeV)", "Events", rebin=bins)
    
    #singlePlot({"name": "recoil_corr_perp", "axis": "recoil_perp" }, "recoil_corr_perp_nolog", -200, 200, 0, 5000, "U_{#perp}  (GeV)", "Events", rebin=5, logY=False)
    
    
    
    