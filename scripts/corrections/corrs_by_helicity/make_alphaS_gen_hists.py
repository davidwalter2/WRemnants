"""
Convenience script to generate the alphaS gen histograms fora given set of theory corrections.
"""

import argparse
import os
from datetime import datetime

THEORY_PREDS = {
    "scetlib_dyturbo_CT18Z_N3p0LL_N2LO_pdfas": {"pdf": "ct18z"},
    "scetlib_dyturbo_CT18Z_N3p1LL_N2LO_pdfas": {"pdf": "ct18z"},
    "scetlib_dyturbo_CT18Z_N4p0LL_N2LO_pdfas": {"pdf": "ct18z"},
    "scetlib_dyturbo_MSHT20_N3p0LL_N2LO_pdfas": {"pdf": "msht20"},
    "scetlib_dyturbo_MSHT20an3lo_N3p0LL_N2LO_pdfas": {"pdf": "msht20an3lo"},
    "scetlib_dyturbo_LatticeNP_CT18Z_N2p1LL_N2LO_pdfas": {"pdf": "ct18z"},
    "scetlib_dyturbo_LatticeNP_CT18Z_N3p0LL_N2LO_pdfas": {"pdf": "ct18z"},
    "scetlib_dyturbo_LatticeNP_CT18Z_N3p1LL_N2LO_pdfas": {"pdf": "ct18z"},
    "scetlib_dyturbo_LatticeNP_CT18Z_N4p0LL_N2LO_pdfas": {"pdf": "ct18z"},
    "scetlib_nnlojet_CT18Z_N3p1LL_N3LO_pdfas": {"pdf": "ct18z"},
    "scetlib_nnlojet_CT18Z_N4p0LL_N3LO_pdfas": {"pdf": "ct18z"},
    "scetlib_nnlojet_MSHT20an3lo_N4p0LL_N3LO_pdfas": {"pdf": "msht20an3lo"},
    "scetlib_dyturbo_LatticeNP_CT18Z_N2p1LL_N2L0_pdfas": {"pdf": "ct18z"},
    "scetlib_dyturbo_LatticeNP_CT18Z_N3p1LL_N2L0_pdfas": {"pdf": "ct18z"},
    "scetlib_dyturbo_LatticeNP_CT18Z_N4p0LL_N2L0_pdfas": {"pdf": "ct18z"}
}


def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--preds",
        nargs="+",
        help="List of theory preds to process. (default: %(default)s)",
        default=list(THEORY_PREDS.keys()),
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        default=f"{os.environ.get('MY_OUT_DIR', '.')}/{datetime.now().strftime('%y%m%d')}_gen_alphaSByHelicity/",
    )
    parser.add_argument(
        "--skim",
        action="store_true",
        help="If set, will run a skimming step to only keep the PDF histograms in the file, saving a new output file.",
    )
    parser.add_argument(
        "--bosons",
        nargs="+",
        choices=["W", "Z"],
        default=["W", "Z"],
        help="Bosons to process. Choose one or both of: W Z (default: %(default)s)",
    )

    return parser.parse_args()


def main():

    args = parse_arguments()

    print("Generating histograms by helicity for the following theory preds:")
    print(args.preds)
    print("Processing bosons:")
    print(args.bosons)
    print("Will output to directory:")
    print(args.outdir)

    filter_procs = []
    aggregate_groups = []
    if "Z" in args.bosons:
        filter_procs.append("Zmumu_13TeVGen")
        aggregate_groups.append("Zmumu")
    if "W" in args.bosons:
        filter_procs.extend(["Wplusmunu_13TeVGen", "Wminusmunu_13TeVGen"])
        aggregate_groups.append("Wmunu")

    filter_procs_arg = " ".join(f"'{proc}'" for proc in filter_procs)
    aggregate_groups_arg = " ".join(aggregate_groups)

    for pred in args.preds:

        pdf = THEORY_PREDS[pred]["pdf"]

        command = f"""
            python {os.environ['WREM_BASE']}/scripts/histmakers/w_z_gen_dists.py --theoryCorr {pred} \
            --filterProcs {filter_procs_arg} --aggregateGroups {aggregate_groups_arg} \
            -o {args.outdir} --maxFiles -1 -j 300 --addHelicityAxis --pdf {pdf}
            """
        print(f"Running command: {command}")
        os.system(command)

        if args.skim:
            pdf_replace = f"_{pdf}" if pdf != "ct18z" else ""
            skim_command = f"""
            python {os.environ['WREM_BASE']}/utilities/open_narf_h5py.py {args.outdir}/w_z_gen_dists_{pred + "_Corr"}_maxFiles_m1{pdf_replace}.hdf5 \
            --filterHistsRegex '^(.*pdfas.*|nominal_gen_theory_uncorr)$' --outfile {args.outdir}/w_z_gen_dists_{pred + "_Corr"}_maxFiles_m1_skimmed.hdf5
            """
            print(f"Running skimming command: {skim_command}")
            os.system(skim_command)

    print("All done!")


if __name__ == "__main__":
    main()
