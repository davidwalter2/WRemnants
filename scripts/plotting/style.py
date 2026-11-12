import matplotlib as mpl

from wums import boostHistHelpers as hh
from wums import logging

logger = logging.child_logger(__name__)


axis_labels = {
    "time": {"label": r"Sidereal time", "unit": "h"},
    "pt_probe": {"label": r"$p_{\text{T}}$", "unit": "GeV"},
    "mll": {"label": r"$m_{\ell \ell}$", "unit": "GeV"},
    "eta_probe": {"label": r"$\eta$", "unit": ""},
}

process_colors = {
    "stat": "red",
    "linearity": "blue",
    "stability": "orange",
    "prefiring_syst": "purple",
    "prefiring_stat": "grey",
}


nuisance_grouping = {
    "liv_unc": [
        "stat",
        "linearity",
        "stability",
        "prefiring_syst",
        "prefiring_stat",
        "bkg",
    ],
    "backgrounds": [
        "Diboson",
        "GG",
        "QCD",
        "QG_2L",
        "QG_Lnu",
        "Top",
        "W_minus",
        "Ztautau",
    ],
}


process_labels = {
    "Data": "Data",
    "Zmumu": r"$Z \rightarrow \mu \mu$",  # $Z \rightarrow \mu \mu $",
    "Diboson": "Diboson",
    "GG": "$\gamma \gamma$",
    "QCD": "QCD",
    "QG_2L": r"$q \gamma \rightarrow \ell \ell$",
    "QG_Lnu": r"$q \gamma \rightarrow \ell \nu$",
    "Top": "Top",
    "W_minus": r"$W^{-}$",
    "Ztautau": r"$Z \rightarrow \tau\tau$",
}


process_supergroups = {
    "nominal": {
        "Zmumu": [
            "Zmumu",
        ],
        "Backgrounds": [
            "Diboson",
            "GG",
            "QCD",
            "QG_2L",
            "QG_Lnu",
            "Top",
            "W_minus",
            "Ztautau",
        ],
    }
}


def get_labels_colors_procs_sorted(procs):
    # order of the processes in the plots by this list
    procs_sort = [
        "Zmumu",
        "backgroumds",
    ][::-1]

    cmap = mpl.colormaps["tab10"]

    procs = sorted(
        procs, key=lambda x: procs_sort.index(x) if x in procs_sort else len(procs_sort)
    )
    logger.debug(f"Found processes {procs} in fitresult")
    labels = [process_labels.get(p, p) for p in procs]
    colors = [process_colors.get(p, cmap(i % cmap.N)) for i, p in enumerate(procs)]
    return labels, colors, procs


def process_grouping(grouping, hist_stack, procs):
    if grouping in process_supergroups.keys():
        new_stack = {}
        for new_name, old_procs in process_supergroups[grouping].items():
            stacks = [hist_stack[procs.index(p)] for p in old_procs if p in procs]
            if len(stacks) == 0:
                continue
            new_stack[new_name] = hh.sumHists(stacks)
    else:
        new_stack = hist_stack
        logger.warning(
            f"No supergroups found for input file with mode {grouping}, proceed without merging groups"
        )

    labels, colors, procs = get_labels_colors_procs_sorted(
        [k for k in new_stack.keys()]
    )
    hist_stack = [new_stack[p] for p in procs]

    return hist_stack, labels, colors, procs
