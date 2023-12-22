import ROOT
import narf
import pathlib

ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")

narf.clingutils.Declare('#include "definitions.h"')
narf.clingutils.Declare('#include "utilities.h"')
narf.clingutils.Declare('#include "pairing.h"')


def make_collection(df, name, selection, input_name="Muon", subselection=None, do_gen_match=False):

    selection = f"[{selection}]"
    if subselection:
        selection += f"[{subselection}]"

    df = df.Define(f"{name}_pt", f"{input_name}_pt{selection}")
    df = df.Define(f"{name}_eta", f"{input_name}_eta{selection}")
    df = df.Define(f"{name}_phi", f"{input_name}_phi{selection}")
        
    if input_name == "MergedStandAloneMuon":
        # check if reco muon has a gen match; for standalone muons within dR<0.3
        if do_gen_match:
            df = df.Define(f"{name}_isGen", f"""
                (zcount::hasMatchDR2(
                    {name}_eta, {name}_phi, GenPart_eta[postFSRmuons][postFSRmuonIdx],GenPart_phi[postFSRmuons][postFSRmuonIdx], 0.09))
                || (zcount::hasMatchDR2(
                    {name}_eta, {name}_phi, GenPart_eta[postFSRantimuons][postFSRantimuonIdx],GenPart_phi[postFSRantimuons][postFSRantimuonIdx], 0.09))
                """)
    else:
        df = df.Define(f"{name}_charge", f"{input_name}_charge{selection}")
        # check if reco muon has a gen match; for regular muons within dR<0.05
        if do_gen_match:
            df = df.Define(f"{name}_isGen", f"""
                ({name}_charge == -1 && zcount::hasMatchDR2(
                    {name}_eta, {name}_phi, GenPart_eta[postFSRmuons][postFSRmuonIdx],GenPart_phi[postFSRmuons][postFSRmuonIdx], 0.0025))
                || ({name}_charge == 1 && zcount::hasMatchDR2(
                    {name}_eta, {name}_phi, GenPart_eta[postFSRantimuons][postFSRantimuonIdx],GenPart_phi[postFSRantimuons][postFSRantimuonIdx], 0.0025))
                """)

    df = df.Define(f"{name}_acceptance", f"{input_name}_acceptance{selection}") 
    df = df.Define(f"{name}_mom4", f"zcount::makeLorentzVector({name}_pt, {name}_eta, {name}_phi, zcount::muon_mass)")
    return df

def make_pairs(df, new_name, object1_name, object2_name=None, define_acceptance_reco=False, define_acceptance_gen=False, define_os_charge=False, define_mass=False):
    if define_mass:
        if object2_name:
            df = df.Define(f"{new_name}_mass", f"zcount::makePairsMass({object1_name}_mom4, {object2_name}_mom4)")
        else:
            df = df.Define(f"{new_name}_mass", f"zcount::makePairsMass({object1_name}_mom4)")

    if define_os_charge:
        if object2_name:
            df = df.Define(f"{new_name}_os", f"zcount::makePairsOS({object1_name}_charge, {object2_name}_charge)")
        else:
            df = df.Define(f"{new_name}_os", f"zcount::makePairsOS({object1_name}_charge)")

    if define_acceptance_reco:
        if object2_name:
            df = df.Define(f"{new_name}_acceptance_reco", f"zcount::makePairsSum({object1_name}_acceptance, {object2_name}_acceptance)")
        else:
            df = df.Define(f"{new_name}_acceptance_reco", f"zcount::makePairsSum({object1_name}_acceptance)")

    if define_acceptance_gen:
        if object2_name:
            df = df.Define(f"{new_name}_isGen", f"zcount::makePairsGen({object1_name}_isGen, {object2_name}_isGen)")
        else:
            df = df.Define(f"{new_name}_isGen", f"zcount::makePairsGen({object1_name}_isGen)")

        # gen level acceptance on an pair level (-1: not a signal pair; 0: out of acceptance; 1: endcap-endcap; 2: endcap-barrel; 3: barrel-barrel)
        df = df.Define(f"{new_name}_acceptance_gen", f"((Event_acceptance_gen + 1) * {new_name}_isGen) - 1")

    return df