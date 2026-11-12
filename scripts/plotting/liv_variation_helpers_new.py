import hist
import numpy as np

from wums.boostHistHelpers import (
    addHists,
    broadcastSystHist,
    divideHists,
    multiplyHists,
    scaleHist,
)

mass_bin = 9


def mc_scaling(hist_in, hist_proj, lumi_scaling, weightsum, cross_sec):
    hist_in /= weightsum
    hist_in *= cross_sec
    hist_in *= 1000
    hist_in_2d = broadcastSystHist(hist_in, hist_proj)
    hist_in_2d = multiplyHists(hist_in_2d, lumi_scaling)

    return hist_in_2d


def make_ones_hist(hist_ref):
    ones = np.ones_like(hist_ref.values())
    h_ones = hist_ref.copy()
    h_ones.values()[...] = ones
    return h_ones


##### FAST MC CORRECTIONS #####
# The corrections below reproduce mc_scaling -> era sum -> make_mutually_exclusive
# with plain NumPy instead of chained boost-histogram operations on the
# (time, mll, pt_probe, eta_probe, pt_tag, eta_tag) hists, which dominated the
# run time. The MC hists have no time axis and every correction is a per-era
# time weight, so val(t, x) = sum_era w_era(t) * mc_era(x). The variance of
# multiplyHists also factorises:
#   (v1 v2)^2 (rel1 + rel2) = v2^2 * [v1^2 rel1] + [v2^2 rel2] * v1^2
# which lets the tag axes be summed *before* the time axis is broadcast.

VAR_CUTOFF = 1e-5  # relVariance cutoff used by multiplyHists
MC_AXES = ("mll", "pt_probe", "eta_probe", "pt_tag", "eta_tag")
PROBE_AXES = ("time", "mll", "pt_probe", "eta_probe")


def _is_weighted(h):
    return h.storage_type == hist.storage.Weight


def _rel_variance(vals, variances):
    return variances / np.clip(vals * vals, VAR_CUTOFF, None)


def _inner(ax, start=0):
    """Slice of the non-flow bins of an axis in a flow=True array."""
    lo = int(ax.traits.underflow)
    return slice(lo + start, ax.extent - int(ax.traits.overflow))


def _make_hist(axes, vals, variances=None):
    if variances is None:
        h = hist.Hist(*axes)
    else:
        h = hist.Hist(*axes, storage=hist.storage.Weight())
        h.variances(flow=True)[...] = variances
    h.values(flow=True)[...] = vals
    return h


def _era_time_weights(lumi_hists, scaling):
    """(values, relative variances) of the H and BG time weights of get_mc_lumis."""
    lumi_h, lumi_bg = lumi_hists
    sum_lumis = addHists(lumi_bg, lumi_h)
    weights = []
    for lumi in (lumi_h, lumi_bg):
        w = multiplyHists(divideHists(lumi, sum_lumis), scaling)
        vals = w.values(flow=True)
        rel = _rel_variance(vals, w.variances(flow=True)) if _is_weighted(w) else None
        weights.append((vals, rel))
    return weights


def _scaled_mc(h, weightsum, cross_sec):
    """Values and variances (None if unweighted) of a normalised 5D MC hist,
    with flow and axes ordered as MC_AXES."""
    h = h.copy()
    h /= weightsum
    h *= cross_sec
    h *= 1000
    order = [h.axes.name.index(n) for n in MC_AXES]
    vals = np.transpose(h.values(flow=True), order)
    variances = np.transpose(h.variances(flow=True), order) if _is_weighted(h) else None
    return vals, variances


def _check_time_axes(time_proj_low):
    if tuple(time_proj_low.axes.name) != ("time", *MC_AXES):
        raise ValueError(f"Unexpected axes {time_proj_low.axes.name}")
    if time_proj_low.axes["time"].traits.underflow or (
        time_proj_low.axes["time"].traits.overflow
    ):
        raise ValueError("time axis with flow bins is not supported")


def get_mc_lumis(
    input_data,
    time_proj_low,
    scaling,
    lumi_hists,
    weightsum,
    cross_sec,
):
    """Normalised, luminosity-weighted MC for (iso, dtdt, dtst, stst), summed over eras.

    input_data is (iso_H, dtdt_H, dtst_H, stst_H, iso_BG, dtdt_BG, dtst_BG, stst_BG).
    Returns full (time, mll, pt_probe, eta_probe, pt_tag, eta_tag) hists.
    """
    _check_time_axes(time_proj_low)
    (w_h, rel_h), (w_bg, rel_bg) = _era_time_weights(lumi_hists, scaling)
    expand = (slice(None),) + (None,) * len(MC_AXES)

    outputs = []
    for h_hist, bg_hist in zip(input_data[:4], input_data[4:]):
        # same order as addHists(bg, h) in the original implementation
        total_vals = total_vars = None
        weighted = True
        for mc_hist, w, rel_w in ((bg_hist, w_bg, rel_bg), (h_hist, w_h, rel_h)):
            v1, var1 = _scaled_mc(mc_hist, weightsum, cross_sec)
            vals = v1[None, ...] * w[expand]
            with_var = var1 is not None and rel_w is not None
            weighted &= with_var
            if with_var:
                variances = vals * vals
                variances *= _rel_variance(v1, var1)[None, ...] + rel_w[expand]
            if total_vals is None:
                total_vals, total_vars = vals, (variances if with_var else None)
            else:
                total_vals += vals
                if weighted:
                    total_vars += variances
        outputs.append(
            _make_hist(time_proj_low.axes, total_vals, total_vars if weighted else None)
        )

    return tuple(outputs)


def _tag_sums(arr, pt_tag_ax, eta_tag_ax):
    """Sum over the trailing tag axes for the tag selections that matter:
    all bins, non-flow bins, and bins with pt_tag >= 25 including the tag
    overflow (the bins remove_bins keeps)."""
    inner = arr[..., _inner(pt_tag_ax), _inner(eta_tag_ax)]
    pt_cut = int(pt_tag_ax.traits.underflow) + pt_tag_ax.index(25)
    return {
        "all": arr.sum(axis=(-2, -1)),
        "inner": inner.sum(axis=(-2, -1)),
        "cut": arr[..., pt_cut:, :].sum(axis=(-2, -1)),
    }


def get_probe_mc_lumis(
    input_data,
    time_proj_low,
    scaling,
    lumi_hists,
    weightsum,
    cross_sec,
):
    """Same as get_mc_lumis + make_mutually_exclusive, followed by
    projecting onto (time, mll, pt_probe, eta_probe), without ever building the
    full tag-resolved hists.

    dtdt is returned with remove_bins already applied (which also drops tag
    bins below 25 GeV before projecting); remove_bins on the result is a no-op.
    """
    _check_time_axes(time_proj_low)
    weights = _era_time_weights(lumi_hists, scaling)
    pt_tag_ax = time_proj_low.axes["pt_tag"]
    eta_tag_ax = time_proj_low.axes["eta_tag"]
    expand = (slice(None), None, None, None)

    regions = []
    weighted = True
    for h_hist, bg_hist in zip(input_data[:4], input_data[4:]):
        total = None
        for mc_hist, (w, rel_w) in ((bg_hist, weights[1]), (h_hist, weights[0])):
            v1, var1 = _scaled_mc(mc_hist, weightsum, cross_sec)
            with_var = var1 is not None and rel_w is not None
            weighted &= with_var
            s_val = _tag_sums(v1, pt_tag_ax, eta_tag_ax)
            if with_var:
                v1sq = v1 * v1
                s_a = _tag_sums(v1sq * _rel_variance(v1, var1), pt_tag_ax, eta_tag_ax)
                s_b = _tag_sums(v1sq, pt_tag_ax, eta_tag_ax)
                c = (w * w * rel_w)[expand]
            era = {}
            for key in s_val:
                vals = w[expand] * s_val[key][None, ...]
                variances = (
                    (w * w)[expand] * s_a[key][None, ...] + c * s_b[key][None, ...]
                    if with_var
                    else None
                )
                era[key] = [vals, variances]
            if total is None:
                total = era
            else:
                for key in total:
                    total[key][0] = total[key][0] + era[key][0]
                    if weighted:
                        total[key][1] = total[key][1] + era[key][1]
        regions.append(total)

    iso, dtdt, dtst, stst = regions
    if not weighted:
        for region in regions:
            for key in region:
                region[key][1] = None

    def _add(a, b, sign):
        if a is None or b is None:
            return None
        return a + b if sign > 0 else a - b

    iso_out = iso["all"]
    dtdt_cut = _add(dtdt["cut"][0], iso["cut"][0], -1)
    # make_mutually_exclusive replaces the first pt_probe bin of dtdt by iso
    # (non-flow bins only) before subtracting it from dtst
    injected = dtdt["all"][0].copy()
    mll_ax, pt_ax, eta_ax = (time_proj_low.axes[n] for n in PROBE_AXES[1:])
    sel = (slice(None), _inner(mll_ax), _inner(pt_ax).start, _inner(eta_ax))
    injected[sel] = iso["inner"][0][sel] + (dtdt["all"][0][sel] - dtdt["inner"][0][sel])
    dtst_out = [
        dtst["all"][0] - injected,
        _add(dtst["all"][1], dtdt["all"][1], 1),
    ]
    stst_out = [
        stst["all"][0] - dtst["all"][0],
        _add(stst["all"][1], dtst["all"][1], 1),
    ]

    axes = [time_proj_low.axes[n] for n in PROBE_AXES]
    return (
        _make_hist(axes, *iso_out),
        remove_bins(_make_hist(axes, dtdt_cut)),
        _make_hist(axes, *dtst_out),
        _make_hist(axes, *stst_out),
    )


def rescale_hists(rescale, hists):
    """Multiply the values of each hist by a time-only hist (e.g. luminometer ratio)."""
    r = rescale.values(flow=True)[(slice(None),) + (None,) * (hists[0].ndim - 1)]
    return tuple(_make_hist(h.axes, r * h.values(flow=True)) for h in hists)


### i need to get good at coding so i dont need to pass in all these variables
def eta_phi_systematic(
    writer,
    input_data,
    time_hists,
    lumi_scaling,
    lumi_hists,
    weightsum,
    cross_sec,
    etaphi_num,
    mass_bin,
):
    (
        iso_H_stat,
        dtdt_H_stat,
        dtst_H_stat,
        stst_H_stat,
        iso_BG_stat,
        dtdt_BG_stat,
        dtst_BG_stat,
        stst_BG_stat,
    ) = input_data

    input_data = [
        iso_H_stat[{"etaPhiRegion": etaphi_num}],
        dtdt_H_stat[{"etaPhiRegion": etaphi_num}],
        dtst_H_stat[{"etaPhiRegion": etaphi_num}],
        stst_H_stat[{"etaPhiRegion": etaphi_num}],
        iso_BG_stat[{"etaPhiRegion": etaphi_num}],
        dtdt_BG_stat[{"etaPhiRegion": etaphi_num}],
        dtst_BG_stat[{"etaPhiRegion": etaphi_num}],
        stst_BG_stat[{"etaPhiRegion": etaphi_num}],
    ]

    iso_stat, dtdt_stat, dtst_stat, stst_stat = get_probe_mc_lumis(
        input_data,
        time_hists,
        lumi_scaling,
        lumi_hists,
        weightsum,
        cross_sec,
    )

    poi_prefire_stat_proj = iso_stat.project("time", "mll")
    writer.add_systematic(
        remove_bins(poi_prefire_stat_proj, "mll", mass_bin + 1),
        f"prefiring_stat_etaphi_{etaphi_num}",
        "Zmumu",
        "ch_iso_poi_high",
        constrained=True,
        groups=["prefiring_stat"],
    )
    writer.add_systematic(
        remove_bins(poi_prefire_stat_proj, "mll", mass_bin + 1, False),
        f"prefiring_stat_etaphi_{etaphi_num}",
        "Zmumu",
        "ch_iso_poi_low",
        constrained=True,
        groups=["prefiring_stat"],
    )

    writer.add_systematic(
        iso_stat[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"prefiring_stat_etaphi_{etaphi_num}",
        "Zmumu",
        "ch_iso_eff",
        constrained=True,
        groups=["prefiring_stat"],
    )
    writer.add_systematic(
        dtdt_stat[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"prefiring_stat_etaphi_{etaphi_num}",
        "Zmumu",
        "ch_dtdt_eff",
        constrained=True,
        groups=["prefiring_stat"],
    )
    writer.add_systematic(
        dtst_stat[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"prefiring_stat_etaphi_{etaphi_num}",
        "Zmumu",
        "ch_dtst_eff",
        constrained=True,
        groups=["prefiring_stat"],
    )
    writer.add_systematic(
        stst_stat[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"prefiring_stat_etaphi_{etaphi_num}",
        "Zmumu",
        "ch_stst_eff",
        constrained=True,
        groups=["prefiring_stat"],
    )


def get_era_vals(mc, trigger_cut, era, type_gen="pass", iso=False):
    if type_gen == "pass":
        if not iso:
            return (
                mc[f"{trigger_cut}_prpg_{era}"].get(),
                mc[f"{trigger_cut}_prpg_{era}_muonL1PrefireSyst"].get(),
                mc[f"{trigger_cut}_prpg_{era}_muonL1PrefireStat"].get(),
            )
        else:
            return (
                mc[f"{trigger_cut}_{era}"].get(),
                mc[f"{trigger_cut}_{era}_muonL1PrefireSyst"].get(),
                mc[f"{trigger_cut}_{era}_muonL1PrefireStat"].get(),
            )

    else:
        return (
            mc[f"{trigger_cut}_prfg_{era}"].get(),
            mc[f"{trigger_cut}_prfg_{era}_muonL1PrefireSyst"].get(),
            mc[f"{trigger_cut}_prfg_{era}_muonL1PrefireStat"].get(),
        )


def luminometer_syst(writer, luminometer, iso, dtdt, dtst, stst, syst):
    dtdt = remove_bins(dtdt)
    writer.add_systematic(
        iso[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"{luminometer}_{syst}",
        "Zmumu",
        "ch_iso_eff",
        constrained=True,
        groups=[f"{syst}"],
    )

    writer.add_systematic(
        dtdt[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"{luminometer}_{syst}",
        "Zmumu",
        "ch_dtdt_eff",
        constrained=True,
        groups=[f"{syst}"],
    )
    writer.add_systematic(
        dtst[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"{luminometer}_{syst}",
        "Zmumu",
        "ch_dtst_eff",
        constrained=True,
        groups=[f"{syst}"],
    )
    writer.add_systematic(
        stst[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"{luminometer}_{syst}",
        "Zmumu",
        "ch_stst_eff",
        constrained=True,
        groups=[f"{syst}"],
    )

    poi_iso_lumi_proj = iso.project("time", "mll")
    writer.add_systematic(
        remove_bins(poi_iso_lumi_proj, "mll", mass_bin + 1),
        f"{luminometer}_{syst}",
        "Zmumu",
        "ch_iso_poi_high",
        constrained=True,
        groups=[f"{syst}"],
    )
    writer.add_systematic(
        remove_bins(poi_iso_lumi_proj, "mll", mass_bin + 1, False),
        f"{luminometer}_{syst}",
        "Zmumu",
        "ch_iso_poi_low",
        constrained=True,
        groups=[f"{syst}"],
    )


def background_syst(
    writer,
    results,
    res_str,
    time_proj_low,
    lumi_scaling,
    lumi_hists,
    proc_name,
    bkg_name,
    fail_gen=False,
    mass_bin=mass_bin,
):

    MC = results[res_str]["output"]
    try:
        weightsum = results[res_str]["weight_sum"]
        cross_sec = results[res_str]["dataset"]["xsec"]
    except:
        weightsum = results["SingleMuon_2016PostVFP"]["weight_sum"]
        cross_sec = results["SingleMuon_2016PostVFP"]["dataset"]["xsec"]
    print(res_str)

    if res_str == "Ztautau_2016PostVFP":
        signal_proc = True
    else:
        signal_proc = False
    ### MAKE THIS IMPLEMENTATION NOT STUPID
    if fail_gen:
        dtdt_prpg_BG, _, _ = get_era_vals(MC, "dtdt", "BG", "fail")
        dtst_prpg_BG, _, _ = get_era_vals(MC, "dtst", "BG", "fail")
        stst_prpg_BG, _, _ = get_era_vals(MC, "stst", "BG", "fail")

        dtdt_prpg_H, _, _ = get_era_vals(MC, "dtdt", "H", "fail")
        dtst_prpg_H, _, _ = get_era_vals(MC, "dtst", "H", "fail")
        stst_prpg_H, _, _ = get_era_vals(MC, "stst", "H", "fail")
    else:
        iso_BG, _, _ = get_era_vals(MC, "pass_iso", "BG", iso=True)
        dtdt_prpg_BG, _, _ = get_era_vals(MC, "dtdt", "BG")
        dtst_prpg_BG, _, _ = get_era_vals(MC, "dtst", "BG")
        stst_prpg_BG, _, _ = get_era_vals(MC, "stst", "BG")

        iso_H, _, _ = get_era_vals(MC, "pass_iso", "H", iso=True)
        dtdt_prpg_H, _, _ = get_era_vals(MC, "dtdt", "H")
        dtst_prpg_H, _, _ = get_era_vals(MC, "dtst", "H")
        stst_prpg_H, _, _ = get_era_vals(MC, "stst", "H")

    prpg_all = [
        iso_H,
        dtdt_prpg_H,
        dtst_prpg_H,
        stst_prpg_H,
        iso_BG,
        dtdt_prpg_BG,
        dtst_prpg_BG,
        stst_prpg_BG,
    ]

    iso, dtdt, dtst, stst = get_probe_mc_lumis(
        prpg_all,
        time_proj_low,
        lumi_scaling,
        lumi_hists,
        weightsum,
        cross_sec,
    )

    iso_poi = iso.project("time", "mll")

    iso_poi_high = remove_bins(iso_poi, "mll", mass_bin + 1)
    writer.add_process(
        iso_poi_high,
        f"{proc_name}",
        "ch_iso_poi_high",
        signal=signal_proc,
    )
    iso_poi_low = remove_bins(
        iso_poi, "mll", mass_bin + 1, False
    )  ### is this the right one?
    writer.add_process(
        iso_poi_low,
        f"{proc_name}",
        "ch_iso_poi_low",
        signal=signal_proc,
    )

    # iso_low = expand_hist by duplicate axes
    # for i in range(nbins_mll):
    #     if i != mass_bin:
    #         if i < mass_bin:
    #             writer.add_systematic(
    #                 addHists(poi_mll_low[{"gen_mll": i}] * var_size, h3_poi_low),
    #                 f"n_mll{i}_{proc_name}",
    #                 f"{proc_name}",
    #                 "ch_iso_poi_low",
    #                 constrained=False,
    #                 groups=["nz"],
    #             )
    #         else:
    #             print("above mass bin")
    #             writer.add_systematic(
    #                 addHists(
    #                     poi_mll_high[{"gen_mll": i - (mass_bin + 1)}] * var_size,
    #                     h3_poi_high,
    #                 ),
    #                 f"n_mll{i}_{proc_name}",
    #                 f"{proc_name}",
    #                 "ch_iso_poi_high",
    #                 constrained=False,
    #                 groups=["nz"],
    #             )

    var = 1.01
    if proc_name in ("W", "W_plus", "W_minus", "Diboson"):
        var = 1.001
    # writer.add_norm_systematic(
    #     f"{bkg_name}", f"{proc_name}", "ch_iso_poi_high", var, groups=["bkg"]
    # )

    # writer.add_norm_systematic(
    #     f"{bkg_name}", f"{proc_name}", "ch_iso_poi_low", var, groups=["bkg"]
    # )

    ##### in efficiency channels ####
    iso_proc = iso[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe")
    dtdt_proc = dtdt[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe")
    dtst_proc = dtst[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe")
    stst_proc = stst[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe")

    if proc_name == "W_minus":
        dtdt_proc.values()[...] = np.abs(
            dtdt_proc.values()
        )  ## there are no variances less than -0.01 so setting this as positive per kenneth's recommendation

    writer.add_process(
        iso_proc, f"n_mll9_{proc_name}", "ch_iso_eff", signal=signal_proc
    )
    writer.add_process(
        dtdt_proc, f"n_mll9_{proc_name}", "ch_dtdt_eff", signal=signal_proc
    )
    writer.add_process(
        dtst_proc, f"n_mll9_{proc_name}", "ch_dtst_eff", signal=signal_proc
    )
    writer.add_process(
        stst_proc, f"n_mll9_{proc_name}", "ch_stst_eff", signal=signal_proc
    )

    var = 1.01
    if proc_name in ("W", "W_plus", "W_minus", "Diboson"):
        var = 1.001

    writer.add_norm_systematic(
        f"{bkg_name}", f"n_mll9_{proc_name}", "ch_iso_eff", var, groups=["bkg"]
    )
    writer.add_norm_systematic(
        f"{bkg_name}", f"n_mll9_{proc_name}", "ch_dtdt_eff", var, groups=["bkg"]
    )
    writer.add_norm_systematic(
        f"{bkg_name}", f"n_mll9_{proc_name}", "ch_dtst_eff", var, groups=["bkg"]
    )
    writer.add_norm_systematic(
        f"{bkg_name}", f"n_mll9_{proc_name}", "ch_stst_eff", var, groups=["bkg"]
    )


def remove_bins(old_hist, ax_name="pt_probe", nbins=1, low=True):
    if ax_name == "pt_probe":
        nbins = old_hist.axes[ax_name].index(25)

    if low:
        new_edges = old_hist.axes[ax_name].edges[nbins:]
    elif not low:
        new_edges = old_hist.axes[ax_name].edges[:nbins]
    new_axis = hist.axis.Variable(new_edges, name=ax_name)
    ax_name_ind = old_hist.axes.name.index(ax_name)

    axes = list(old_hist.axes)
    axes[ax_name_ind] = new_axis

    # slices in flow=True coordinates; flow bins are dropped except on the tag
    # axes, where they are kept (the pt_tag overflow is populated)
    src = [_inner(ax) for ax in old_hist.axes]
    first = src[ax_name_ind].start
    if low:
        src[ax_name_ind] = slice(first + nbins, src[ax_name_ind].stop)
    elif not low:
        src[ax_name_ind] = slice(first, first + nbins - 1)
    dst = [_inner(ax) for ax in axes]
    dst[ax_name_ind] = _inner(new_axis)

    tag_name = f"{ax_name[:-5]}tag"
    if "probe" in ax_name and tag_name in old_hist.axes.name:
        tag_ind = old_hist.axes.name.index(tag_name)
        old_tag_axis = old_hist.axes[tag_ind]
        new_tag_axis = hist.axis.Variable(
            new_edges, name=tag_name, overflow=old_tag_axis.traits.overflow
        )
        axes[tag_ind] = new_tag_axis
        src[tag_ind] = slice(int(old_tag_axis.traits.underflow) + nbins, None)
        dst[tag_ind] = slice(int(new_tag_axis.traits.underflow), None)
        other_tag = "eta_tag" if tag_name == "pt_tag" else "pt_tag"
        if other_tag in old_hist.axes.name:
            other_ind = old_hist.axes.name.index(other_tag)
            src[other_ind] = dst[other_ind] = slice(None)

    new_hist = hist.Hist(*axes)
    new_hist.values(flow=True)[tuple(dst)] = old_hist.values(flow=True)[tuple(src)]

    return new_hist


def make_mutually_exclusive(iso, dtdt, dtst, stst):
    def _subtract(h1, vals2, h2):
        # addHists(h1, scaleHist(h2, -1)), with the values of h2 replaced by vals2
        weighted = _is_weighted(h1) and _is_weighted(h2)
        return _make_hist(
            h1.axes,
            h1.values(flow=True) - vals2,
            h1.variances(flow=True) + h2.variances(flow=True) if weighted else None,
        )

    iso_ex = iso
    dtdt_ex = _subtract(dtdt, iso.values(flow=True), iso)
    # the first pt_probe bin of dtdt is replaced by iso (non-flow bins only)
    dtdt_iso_injected = dtdt.values(flow=True).copy()
    sel = tuple(
        _inner(ax).start if ax.name == "pt_probe" else _inner(ax) for ax in dtdt.axes
    )
    dtdt_iso_injected[sel] = iso.values(flow=True)[sel]
    dtst_ex = _subtract(dtst, dtdt_iso_injected, dtdt)
    stst_ex = _subtract(stst, dtst.values(flow=True), dtst)

    return iso_ex, dtdt_ex, dtst_ex, stst_ex


def create_variation(
    variation_hist,
    refererence_hist,
    i,
    j,
    k,
    nbins_total,
    muon="tag",
    h2=False,
    var_size=0.01,
):

    not_muon = "probe"
    var = variation_hist[{"gen_time": k, f"pt_{not_muon}": i, f"eta_{not_muon}": j}]

    var = var.project("time", "mll", f"pt_{muon}", f"eta_{muon}")
    var = broadcastSystHist(var, refererence_hist)
    var = multiplyHists(scaleHist(var, var_size / (nbins_total)), refererence_hist)
    var = var.project("time", "mll", "pt_probe", "eta_probe")
    return var


def variation_array(values, axes, selections, output_axes):
    """Select named bins and return a NumPy array in a requested axis order."""
    remaining_axes = list(axes.name)
    for axis_name, index in selections.items():
        axis_index = remaining_axes.index(axis_name)
        values = np.take(values, index, axis=axis_index)
        remaining_axes.pop(axis_index)
    return np.transpose(
        values, [remaining_axes.index(axis_name) for axis_name in output_axes]
    )


def variation_values(variation_hist):
    """Return the compact 4D NumPy variation array without duplicate axes."""
    return np.asarray(variation_hist.values())


def create_variation_array(
    variation_hist,
    reference_values,
    i,
    j,
    nbins_total,
    h2=False,
    var_size=0.01,
):
    """Create all generator-time variations as a NumPy array.

    The returned shape is ``(gen_time, time, mll, pt_probe, eta_probe)``.
    ``reference_values`` must already have the latter four axes in that order.
    """
    hist_i = i - 1 if h2 else i
    factors = variation_array(
        variation_hist.values(),
        variation_hist.axes,
        {"pt_probe": hist_i, "eta_probe": j},
        ["gen_time", "time", "mll", "pt_tag", "eta_tag"],
    )
    if factors.shape[1:] != reference_values.shape:
        raise ValueError(
            f"Variation shape {factors.shape[1:]} does not match reference shape "
            f"{reference_values.shape}"
        )
    return factors * reference_values[None, ...] * (var_size / nbins_total)


def create_probe_variation_array(variation_hist, reference_values, i, j, var_size=0.01):
    """Create all generator-time probe variations as a NumPy array."""
    factors = variation_array(
        variation_hist.values(),
        variation_hist.axes,
        {},
        ["gen_time", "time", "mll", "pt_probe", "eta_probe", "pt_tag", "eta_tag"],
    )[:, :, :, :, :, i, j]
    return factors * reference_values[None, ...] * var_size


def prefiring_syst(writer, iso_prefire, dtdt_prefire, dtst_prefire, stst_prefire):
    writer.add_systematic(
        iso_prefire[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"prefiring_syst",
        "Zmumu",
        "ch_iso_eff",
        constrained=True,
        groups=["prefiring_syst"],
    )
    dtdt_prefire = remove_bins(dtdt_prefire)
    writer.add_systematic(
        dtdt_prefire[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"prefiring_syst",
        "Zmumu",
        "ch_dtdt_eff",
        constrained=True,
        groups=["prefiring_syst"],
    )
    writer.add_systematic(
        dtst_prefire[{"mll": mass_bin}].project("time", "pt_probe", "eta_probe"),
        f"prefiring_syst",
        "Zmumu",
        "ch_dtst_eff",
        constrained=True,
        groups=["prefiring_syst"],
    )
    writer.add_systematic(
        stst_prefire[{"mll": mass_bin}].project(
            "time", "pt_probe", "eta_probe"
        ),  # used to be tag pt and eta
        f"prefiring_syst",
        "Zmumu",
        "ch_stst_eff",
        constrained=True,
        groups=["prefiring_syst"],
    )

    poi_prefire_proj = iso_prefire.project("time", "mll")
    writer.add_systematic(
        remove_bins(poi_prefire_proj, "mll", mass_bin + 1),
        f"prefiring_syst",
        "Zmumu",
        "ch_iso_poi_high",
        constrained=True,
        groups=["prefiring_syst"],
    )

    writer.add_systematic(
        remove_bins(poi_prefire_proj, "mll", mass_bin + 1, False),
        f"prefiring_syst",
        "Zmumu",
        "ch_iso_poi_low",
        constrained=True,
        groups=["prefiring_syst"],
    )


def get_corrected_mc(
    results,
    process,
    time_proj_low,
    lumi_hists,
    pcc_scaling,
    hfoc_scaling,
    ramses_scaling,
    hfoc_sbil,
    ramses_sbil,
    lumi_scaling,
    luminometers=False,
):
    print(process)
    MC = results[process]["output"]

    iso_BG, iso_BG_syst, iso_BG_stat = get_era_vals(MC, "pass_iso", "BG", iso=True)
    dtdt_BG, dtdt_BG_syst, dtdt_BG_stat = get_era_vals(MC, "dtdt", "BG")
    dtst_BG, dtst_BG_syst, dtst_BG_stat = get_era_vals(MC, "dtst", "BG")
    stst_BG, stst_BG_syst, stst_BG_stat = get_era_vals(MC, "stst", "BG")

    iso_H, iso_H_syst, iso_H_stat = get_era_vals(MC, "pass_iso", "H", iso=True)
    dtdt_H, dtdt_H_syst, dtdt_H_stat = get_era_vals(MC, "dtdt", "H")
    dtst_H, dtst_H_syst, dtst_H_stat = get_era_vals(MC, "dtst", "H")
    stst_H, stst_H_syst, stst_H_stat = get_era_vals(MC, "stst", "H")
    syst = [
        iso_H_syst[{"downUpVar": 0}],
        dtdt_H_syst[{"downUpVar": 0}],
        dtst_H_syst[{"downUpVar": 0}],
        stst_H_syst[{"downUpVar": 0}],
        iso_BG_syst[{"downUpVar": 0}],
        dtdt_BG_syst[{"downUpVar": 0}],
        dtst_BG_syst[{"downUpVar": 0}],
        stst_BG_syst[{"downUpVar": 0}],
    ]
    stat = [
        iso_H_stat[{"downUpVar": 0}],
        dtdt_H_stat[{"downUpVar": 0}],
        dtst_H_stat[{"downUpVar": 0}],
        stst_H_stat[{"downUpVar": 0}],
        iso_BG_stat[{"downUpVar": 0}],
        dtdt_BG_stat[{"downUpVar": 0}],
        dtst_BG_stat[{"downUpVar": 0}],
        stst_BG_stat[{"downUpVar": 0}],
    ]

    prpg_all = [
        iso_H,
        dtdt_H,
        dtst_H,
        stst_H,
        iso_BG,
        dtdt_BG,
        dtst_BG,
        stst_BG,
    ]

    weightsum = results[process]["weight_sum"]
    xsec = results[process]["dataset"]["xsec"]

    pass_gen = MC["pass_gen"].get()

    lumi_args = (time_proj_low, lumi_scaling, lumi_hists, weightsum, xsec)

    iso, dtdt, dtst, stst = make_mutually_exclusive(*get_mc_lumis(prpg_all, *lumi_args))
    corrected_mc = [iso, dtdt, dtst, stst]
    print("through nominal corrections")

    # the variations below are only used through projections onto the probe
    # axes, so they are built directly as (time, mll, pt_probe, eta_probe) hists
    if luminometers:
        probe_mc = get_probe_mc_lumis(prpg_all, *lumi_args)

        hfoc_stability = rescale_hists(
            divideHists(hfoc_scaling, lumi_scaling), probe_mc
        )
        pcc_stability = rescale_hists(divideHists(pcc_scaling, lumi_scaling), probe_mc)
        ramses_stability = rescale_hists(
            divideHists(ramses_scaling, lumi_scaling), probe_mc
        )
        hfoc_linearity = rescale_hists(divideHists(hfoc_sbil, lumi_scaling), probe_mc)
        ramses_linearity = rescale_hists(
            divideHists(ramses_sbil, lumi_scaling), probe_mc
        )
        print("luminometers")

    ### dtdt was negative from before this was all passed into a single function. need to investigate further
    iso_prefire, dtdt_prefire, dtst_prefire, stst_prefire = get_probe_mc_lumis(
        syst, *lumi_args
    )
    print("prefiring")

    pass_gen = mc_scaling(
        pass_gen,
        time_proj_low,
        lumi_scaling,
        weightsum,
        xsec,
    )
    corrected_prefiring = [iso_prefire, dtdt_prefire, dtst_prefire, stst_prefire]

    if luminometers:
        return (
            corrected_mc,
            pass_gen,
            corrected_prefiring,
            stat,
            weightsum,
            xsec,
            pcc_stability,
            hfoc_stability,
            ramses_stability,
            hfoc_linearity,
            ramses_linearity,
        )

    else:
        return corrected_mc, pass_gen, corrected_prefiring, stat
