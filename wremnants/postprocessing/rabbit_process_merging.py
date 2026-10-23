"""Write groups of processes to the rabbit tensor as one process each."""

import hist
import numpy as np

from wums import logging

from wums import boostHistHelpers as hh  # isort: skip

logger = logging.child_logger(__name__)

# W = Wmunu + Wtaunu, Z = Zmumu + Ztautau, TopVV = Top + Diboson
default_merges = {
    "W": ["Wmunu", "Wtaunu"],
    "Z": ["Zmumu", "Ztautau"],
    "TopVV": ["Top", "Diboson"],
}


def parse_merges(specs):
    """{merged process: [components]} from 'NAME=PROC1,PROC2' strings; the default
    merges without any."""
    if not specs:
        return {k: list(v) for k, v in default_merges.items()}
    merges = {}
    for spec in specs:
        name, sep, components = spec.partition("=")
        components = [c for c in components.split(",") if c]
        if not sep or not name or len(components) < 2:
            raise ValueError(
                f"--mergeProcs entries must be of the form NAME=PROC1,PROC2[,...], got '{spec}'"
            )
        merges[name] = components
    return merges


class ProcessMergingWriter:
    """Wraps a rabbit TensorWriter and writes groups of processes as one process
    each, e.g. W = Wmunu + Wtaunu, so that the MC statistics of the components are
    pooled into one template per bin: fewer and larger templates for the
    bin-by-bin statistical uncertainties and the MC-stat de-biasing.

    The nominal templates are summed, with their MC-stat folds. A systematic keeps
    acting on the components it was added for: the merged variation is the sum of
    their variations and of the nominal templates of the other components, so the
    normalisation uncertainty of one component becomes a shape variation of the
    merged process. Everything else is passed through to the wrapped writer.

    The components of a merged process are collected per channel and written by
    finalize_merged_processes(channel), which Datagroups.addNominalHistograms
    calls after adding its processes. The variations of a systematic are
    collected and written by flush_merged_systematics(), which
    Datagroups.addSystematic calls at its end: all components that a systematic
    varies must be varied in the same call."""

    def __init__(self, writer, merges):
        self.writer = writer
        self.merges = {name: list(components) for name, components in merges.items()}
        self.merged_of = {}
        for name, components in self.merges.items():
            for c in components:
                if c in self.merged_of:
                    raise ValueError(
                        f"Process {c} is a component of both {self.merged_of[c]} and {name}"
                    )
                if c in self.merges:
                    raise ValueError(f"Process {c} is both merged and a component")
                self.merged_of[c] = name
        # (channel, merged) -> {component: (hist, fold axis, signal)}, in arrival order
        self.components = {}
        # (channel, merged) -> the fold axis it was written with (None: without folds)
        self.finalized = {}
        # (channel, merged, systematic) -> (kwargs, {component: (hists, fold axis)})
        self.pending = {}
        self.flushed = set()

    def __getattr__(self, name):
        return getattr(self.writer, name)

    def _check_kwargs(self, kwargs):
        if kwargs.get("as_difference", False):
            raise NotImplementedError(
                "Systematics given as differences are not supported for merged processes"
            )
        if kwargs.get("syst_axes"):
            raise NotImplementedError(
                "Explicit systematic axes are not supported for merged processes"
            )

    def add_process(
        self, h, name, channel="ch0", signal=False, variances=None, fold_axis=None
    ):
        if name not in self.merged_of:
            return self.writer.add_process(
                h,
                name,
                channel,
                signal=signal,
                variances=variances,
                fold_axis=fold_axis,
            )
        if variances is not None:
            raise NotImplementedError(
                f"Explicit variances are not supported for the component {name} of a merged process"
            )
        merged = self.merged_of[name]
        if (channel, merged) in self.finalized:
            raise RuntimeError(
                f"Process {merged} in channel {channel} was written already, its component {name} comes too late"
            )
        self.components.setdefault((channel, merged), {})[name] = (h, fold_axis, signal)

    def finalize_merged_processes(self, channel):
        for (ch, merged), comps in self.components.items():
            if ch != channel or (ch, merged) in self.finalized:
                continue
            signals = {s for _, _, s in comps.values()}
            if len(signals) > 1:
                raise RuntimeError(
                    f"The components {list(comps)} of {merged} are not all signal or all background"
                )
            hists = [h for h, _, _ in comps.values()]
            fold_axes = [fa for _, fa, _ in comps.values()]
            fold_axis = fold_axes[0] if len(set(fold_axes)) == 1 else None
            if fold_axis is None and any(fa is not None for fa in fold_axes):
                logger.warning(
                    f"Not all components of {merged} have MC-stat folds, it is written without"
                )
            h = _sum_hists(hists, fold_axes, fold_axis)
            logger.info(
                f"Add process {merged} in channel {channel}, merged from {list(comps)}"
            )
            self.writer.add_process(
                h, merged, channel, signal=signals.pop(), fold_axis=fold_axis
            )
            self.finalized[(ch, merged)] = fold_axis

    def add_systematic(self, h, name, process, channel, **kwargs):
        if process not in self.merged_of:
            return self.writer.add_systematic(h, name, process, channel, **kwargs)
        self._check_kwargs(kwargs)
        merged = self.merged_of[process]
        if (channel, merged) not in self.finalized:
            self.finalize_merged_processes(channel)
        key = (channel, merged, name)
        if key in self.flushed:
            raise RuntimeError(
                f"Systematic {name} of {merged} was written already: the components of a merged process must be varied in the same call"
            )
        fold_axis = kwargs.pop("fold_axis", None)
        stored_kwargs, variations = self.pending.setdefault(key, (kwargs, {}))
        if stored_kwargs != kwargs:
            raise RuntimeError(
                f"Systematic {name} is added with different options for the components of {merged}: {stored_kwargs} and {kwargs}"
            )
        variations[process] = (h, fold_axis)

    def flush_merged_systematics(self):
        for (channel, merged, name), (kwargs, variations) in self.pending.items():
            comps = self.components[(channel, merged)]
            unknown = [c for c in variations if c not in comps]
            if unknown:
                raise RuntimeError(
                    f"Systematic {name} for {unknown}, which are not processes of channel {channel}"
                )
            first = next(iter(variations.values()))[0]
            two_sided = not isinstance(first, hist.Hist)
            # the components without this systematic enter with their nominal
            terms, fold_axes = [], []
            for c, (hnom, fa_nom, _) in comps.items():
                if c in variations:
                    hvar, fa = variations[c]
                    if two_sided == isinstance(hvar, hist.Hist):
                        raise RuntimeError(
                            f"Systematic {name}: the components of {merged} mix one- and two-sided variations"
                        )
                    terms.append(list(hvar) if two_sided else [hvar])
                    fold_axes.append(fa)
                else:
                    terms.append([hnom, hnom] if two_sided else [hnom])
                    fold_axes.append(fa_nom)
            # split per fold only if all terms have the folds, as the merged process
            fold_axis = fold_axes[0] if len(set(fold_axes)) == 1 else None
            if fold_axis is not None and self.finalized[(channel, merged)] is None:
                fold_axis = None
            hists = [
                _sum_hists([t[i] for t in terms], fold_axes, fold_axis)
                for i in range(len(terms[0]))
            ]
            self.writer.add_systematic(
                hists if two_sided else hists[0],
                name,
                merged,
                channel,
                fold_axis=fold_axis,
                **kwargs,
            )
            self.flushed.add((channel, merged, name))
        self.pending = {}

    def add_norm_systematic(self, name, process, channel, uncertainty, **kwargs):
        processes = [process] if isinstance(process, str) else list(process)
        if not any(p in self.merged_of for p in processes):
            return self.writer.add_norm_systematic(
                name, process, channel, uncertainty, **kwargs
            )
        self.finalize_merged_processes(channel)
        uncertainties = (
            list(uncertainty)
            if isinstance(uncertainty, (list, tuple, np.ndarray))
            else [uncertainty]
        )
        if len(uncertainties) == 1:
            uncertainties = uncertainties * len(processes)
        unc_of = dict(zip(processes, uncertainties))
        # a log-normal uncertainty of the merged process only if all its components
        # have the same one
        out = {}
        for p in processes:
            merged = self.merged_of.get(p, p)
            if merged != p and any(
                unc_of.get(c) != unc_of[p]
                for c in self.components.get((channel, merged), {})
            ):
                raise NotImplementedError(
                    f"Normalisation uncertainty {name} on a part of {merged}: add it through Datagroups.addNormSystematic"
                )
            out.setdefault(merged, unc_of[p])
        return self.writer.add_norm_systematic(
            name, list(out), channel, list(out.values()), **kwargs
        )

    def add_beta_variations(self, h, process, *args, **kwargs):
        if process in self.merged_of:
            raise NotImplementedError(
                f"Beta variations are not supported for the component {process} of a merged process"
            )
        return self.writer.add_beta_variations(h, process, *args, **kwargs)

    def write(self, *args, **kwargs):
        self.flush_merged_systematics()
        missing = [k for k in self.components if k not in self.finalized]
        if missing:
            raise RuntimeError(f"Merged processes {missing} were never written")
        return self.writer.write(*args, **kwargs)


def _sum_hists(hists, fold_axes, fold_axis):
    """Sum of the histograms, with the axes of the first one; the MC-stat fold axis
    is summed over in those that have it unless all share `fold_axis`."""
    if fold_axis is None:
        hists = [
            h if fa is None else h.project(*[n for n in h.axes.name if n != fa])
            for h, fa in zip(hists, fold_axes)
        ]
    names = hists[0].axes.name
    out = hists[0].copy()
    for h in hists[1:]:
        if h.axes.name != names:
            if set(h.axes.name) != set(names):
                raise RuntimeError(
                    f"Cannot sum histograms with axes {names} and {h.axes.name}"
                )
            h = h.project(*names)
        out = hh.addHists(out, h, createNew=False)
    return out
