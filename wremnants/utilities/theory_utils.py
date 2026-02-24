import glob
import os
import re

from wremnants.utilities import common

pdfMap = {
    "nnpdf31": {
        "name": "pdfNNPDF31",
        "branch": "LHEPdfWeight",
        "combine": "symHessian",
        "entries": 101,
        "alphas": ["LHEPdfWeight[0]", "LHEPdfWeight[101]", "LHEPdfWeight[102]"],
        "alphasRange": "002",
        "inflation_factor_wmass": 3.0,
        "inflation_factor_alphaS": 3.0,
    },
    "ct18": {
        "name": "pdfCT18",
        "branch": "LHEPdfWeightAltSet11",
        "combine": "asymHessian",
        "entries": 59,
        "alphas": [
            "LHEPdfWeightAltSet11[0]",
            "LHEPdfWeightAltSet11[59]",
            "LHEPdfWeightAltSet11[62]",
        ],
        "alphasRange": "002",
        "scale": 1 / 1.645,  # Convert from 90% CL to 68%
        "inflation_factor_wmass": 1.0,
        "inflation_factor_alphaS": 1.2,
    },
    "nnpdf30": {
        "name": "pdfNNPDF30",
        "branch": "LHEPdfWeightAltSet7",
        "combine": "symHessian",
        "entries": 101,
        "alphas": [
            "LHEPdfWeightAltSet13[0]",
            "LHEPdfWeightAltSet15[0]",
            "LHEPdfWeightAltSet16[0]",
        ],
        "alphasRange": "001",
        "inflation_factor_wmass": 1.0,  # not determined
        "inflation_factor_alphaS": 1.0,  # not determined
    },
    "nnpdf40": {
        "name": "pdfNNPDF40",
        "branch": "LHEPdfWeightAltSet3",
        "combine": "symHessian",
        "entries": 51,
        "alphas": [
            "LHEPdfWeightAltSet3[0]",
            "LHEPdfWeightAltSet3[51]",
            "LHEPdfWeightAltSet3[52]",
        ],
        "alphasRange": "001",
        "inflation_factor_wmass": 5.0,
        "inflation_factor_alphaS": 5.0,
    },
    "pdf4lhc21": {
        "name": "pdfPDF4LHC21",
        "branch": "LHEPdfWeightAltSet10",
        "combine": "symHessian",
        "entries": 41,
        "alphas": [
            "LHEPdfWeightAltSet10[0]",
            "LHEPdfWeightAltSet10[41]",
            "LHEPdfWeightAltSet10[42]",
        ],
        "alphasRange": "001",
        "inflation_factor_wmass": 1.0,
        "inflation_factor_alphaS": 1.5,
    },
    "msht20": {
        "name": "pdfMSHT20",
        "branch": "LHEPdfWeightAltSet12",
        "combine": "asymHessian",
        "entries": 65,
        "alphas": [
            "LHEPdfWeightAltSet12[0]",
            "LHEPdfWeightAltSet12[67]",
            "LHEPdfWeightAltSet12[70]",
        ],
        "alphasRange": "002",
        "inflation_factor_wmass": 1.5,
        "inflation_factor_alphaS": 2.0,
    },
    "msht20mcrange": {
        "name": "pdfMSHT20mcrange",
        "branch": "LHEPdfWeightAltSet12",
        "combine": "asymHessian",
        "entries": 9,
        "first_entry": 72,
    },
    "msht20mbrange": {
        "name": "pdfMSHT20mbrange",
        "branch": "LHEPdfWeightAltSet12",
        "combine": "asymHessian",
        "entries": 7,
        "first_entry": 81,
    },
    "msht20mcrange_renorm": {
        "name": "pdfMSHT20mcrange",
        "branch": "LHEPdfWeightAltSet12",
        "combine": "asymHessian",
        "entries": 9,
        "first_entry": 72,
        "renorm": True,
    },
    "msht20mbrange_renorm": {
        "name": "pdfMSHT20mbrange",
        "branch": "LHEPdfWeightAltSet12",
        "combine": "asymHessian",
        "entries": 7,
        "first_entry": 81,
        "renorm": True,
    },
    "msht20an3lo": {
        "name": "pdfMSHT20an3lo",
        "branch": "LHEPdfWeightAltSet24",
        "combine": "asymHessian",
        "entries": 105,
        "alphas": [
            "LHEPdfWeightAltSet24[0]",
            "LHEPdfWeightAltSet24[108]",
            "LHEPdfWeightAltSet24[111]",
        ],
        "alphasRange": "002",
        "inflation_factor_wmass": 1.5,
        "inflation_factor_alphaS": 1.0,  # not determined
    },
    "ct18z": {
        "name": "pdfCT18Z",
        "branch": "LHEPdfWeightAltSet11",
        "combine": "asymHessian",
        "entries": 59,
        "first_entry": 63,
        "alphas": [
            "LHEPdfWeightAltSet11[63]",
            "LHEPdfWeightAltSet11[122]",
            "LHEPdfWeightAltSet11[125]",
        ],
        "alphasRange": "002",
        "scale": 1 / 1.645,  # Convert from 90% CL to 68%
        "inflation_factor_wmass": 1.0,
        "inflation_factor_alphaS": 1.0,
    },
    "atlasWZj20": {
        "name": "pdfATLASWZJ20",
        "branch": "LHEPdfWeightAltSet19",
        "combine": "asymHessian",
        "entries": 60,
        "alphas": ["LHEPdfWeight[0]", "LHEPdfWeight[41]", "LHEPdfWeight[42]"],
        "alphasRange": "002",
        "inflation_factor_wmass": 1.0,  # not determined
        "inflation_factor_alphaS": 1.0,  # not determined
    },
    "herapdf20": {
        "name": "pdfHERAPDF20",
        "branch": "LHEPdfWeightAltSet20",
        "combine": "asymHessian",
        "entries": 29,
        "alphas": [
            "LHEPdfWeightAltSet20[0]",
            "LHEPdfWeightAltSet22[0]",
            "LHEPdfWeightAltSet23[0]",
        ],  # alphas 116-120
        "alphasRange": "002",
        "inflation_factor_wmass": 4.0,
        "inflation_factor_alphaS": 3.5,
    },
    "herapdf20ext": {
        "name": "pdfHERAPDF20ext",
        "branch": "LHEPdfWeightAltSet21",
        "combine": "asymHessian",
        "entries": 14,
        "alphas": [
            "LHEPdfWeightAltSet20[0]",
            "LHEPdfWeightAltSet22[0]",
            "LHEPdfWeightAltSet23[0]",
        ],  # dummy AS
        "alphasRange": "002",
        "inflation_factor_wmass": 4.0,
        "inflation_factor_alphaS": 3.0,
    },
}


only_central_pdf_datasets = [
    "Wplusmunu_bugfix",
    "Wminusmunu_bugfix",
    "Zmumu_bugfix",
    "Zmumu_bugfix_slc7",
]


def pdf_info_map(dataset, pdfset):
    infoMap = pdfMap

    # Just ignore PDF variations for non W/Z samples
    if (
        pdfset is None
        or not (dataset[0] in ["W", "Z"] and dataset[1] not in ["W", "Z"])
        or "horace" in dataset
        or (pdfset != "nnpdf31" and dataset in only_central_pdf_datasets)
        or pdfset not in infoMap
    ):
        raise ValueError(f"Skipping PDF {pdfset} for dataset {dataset}")
    return infoMap[pdfset]


def pdfNamesSymHessian(entries, pdfset=""):
    return [f"pdf{i+1}{pdfset.replace('pdf', '')}" for i in range(entries)]


def pdfNamesAsymHessian(entries, pdfset=""):
    pdfNames = ["pdf0" + pdfset.replace("pdf", "")]
    if pdfset == "pdfHERAPDF20ext":
        entries -= 3
    pdfNames.extend(
        [
            f"pdf{int((j+2)/2)}{pdfset.replace('pdf', '')}{'Up' if j % 2 else 'Down'}"
            for j in range(entries - 1)
        ]
    )
    if pdfset == "pdfHERAPDF20ext":
        pdfNames.extend(
            [f"pdf{entries//2 + j}{pdfset.replace('pdf', '')}" for j in range(1, 4)]
        )
    return pdfNames


def valid_theory_corrections():
    corr_files = glob.glob(common.data_dir + "TheoryCorrections/*Corr*.pkl.lz4")
    matches = [
        re.match(r"(^.*)_Corr[W|Z|BSM]\.pkl\.lz4", os.path.basename(c))
        for c in corr_files
    ]
    return [m[1] for m in matches if m] + ["none"]


def valid_ew_theory_corrections():
    corr_files = glob.glob(
        common.data_dir + "TheoryCorrections/*[eE][wW]*Corr*.pkl.lz4"
    )
    matches = [
        re.match(r"(^.*)Corr[W|Z]\.pkl\.lz4", os.path.basename(c)) for c in corr_files
    ]
    return [m[1] for m in matches if m] + ["none"]
