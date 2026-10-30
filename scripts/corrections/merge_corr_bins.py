"""Adaptive merging of the bins of a MiNNLO reference where the correction by helicity is ill-defined.

The correction by helicity (CentralCorrByHelicityHelper) weights an event by
sum_h sigma_h^corr P_h(Omega) / sum_h sigma_h^MiNNLO P_h(Omega). In a bin of the MiNNLO reference the denominator is
sigma_UL (1 + cos^2 theta + sum_i A_i P_i(theta, phi)); it is ill-defined where sigma_UL <= 0 or where the bracket is
not positive at all angles, which happens in sparsely populated bins (high |Y|, high pT, high mass).

The reference keeps its rectangular grid; bins are grouped, and the correction is to be computed on the groups
(MiNNLO and the prediction summed over the same groups) and assigned to every bin of a group. A group passes if
sigma_UL > 0, N_eff = sigma_UL^2 / var(sigma_UL) >= nmin and the bracket is >= bmin at all angles. Bins are merged
only where needed, first along the first merge axis (default pT), then along both merge axes (pT and mass); the other
axes (|Y|, charge) are never merged. At each step the worst failing group (sigma_UL <= 0 first, then the smallest
N_eff) is merged with the neighbour that gives a passing group of the smallest N_eff, or else with the neighbour of
the smallest N_eff. Empty bins (no events) are not merged for their own sake; a failing group crosses them
only when it has no populated neighbour. Bins outside the generated mass range of the samples (--massRange) are not
merged at all: MiNNLO is empty there because of the generation cut, and a group extending there would add predicted
cross section without MiNNLO events to its ratio.

Output: an npz file with the group index of every bin (array over the non-helicity axes) and the axis edges.
"""

import argparse
import itertools

import numpy as np

from wremnants.utilities.io_tools import input_tools
from wums import logging

NCOS, NPHI = 61, 96
_c = np.linspace(-1, 1, NCOS)[:, None]
_phi = np.linspace(0, 2 * np.pi, NPHI, endpoint=False)[None, :]
_s = np.sqrt(1 - _c**2)
# angular polynomials of csAngularFactors for A0..A7, on a (cos theta, phi) grid: (8, npoints)
BASIS = np.stack(
    [
        0.5 * (1 - 3 * _c**2) + 0 * _phi,
        2 * _s * _c * np.cos(_phi),
        0.5 * _s**2 * np.cos(2 * _phi),
        _s * np.cos(_phi),
        _c + 0 * _phi,
        _s**2 * np.sin(2 * _phi),
        2 * _s * _c * np.sin(_phi),
        _s * np.sin(_phi),
    ]
).reshape(8, -1)
BASE = ((1 + _c**2) + 0 * _phi).reshape(-1)


def group_stats(s, var):
    """(sigma_UL, N_eff, minimum of the angular bracket) of summed helicity cross sections s and UL variance var"""
    ul = s[0]
    if ul <= 0:
        return ul, 0.0, -np.inf
    neff = ul**2 / var if var > 0 else 0.0
    return ul, neff, float((BASE + (s[1:] / ul) @ BASIS).min())


def passes(x, nmin, bmin):
    return x[0] > 0 and x[1] >= nmin and x[2] >= bmin


def is_empty(s, var):
    # no events: not failing, crossed only to reach populated bins
    return s[0] == 0 and var == 0


def merge_plane(vals, ulvar, active, nmin, bmin):
    """Group the cells of a 2D plane (first merge axis, second merge axis); returns the group label of every cell.
    Inactive cells are left as they are."""
    n0, n1 = vals.shape[:2]
    ncell = n0 * n1
    parent = np.arange(ncell)
    # cells are numbered with the first merge axis running fastest, c = i0 + n0 * i1
    ssum = vals.reshape(ncell, -1, order="F").astype(float)
    svar = ulvar.reshape(-1, order="F").astype(float)
    active = active.reshape(-1, order="F")
    members = {c: [c] for c in range(ncell)}
    st = {c: group_stats(ssum[c], svar[c]) for c in range(ncell)}

    def find(c):
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    def neighbours(root, both_axes):
        out = set()
        steps = (
            ((1, 0), (-1, 0)) if not both_axes else ((1, 0), (-1, 0), (0, 1), (0, -1))
        )
        for c in members[root]:
            i1, i0 = divmod(c, n0)
            for d0, d1 in steps:
                j0, j1 = i0 + d0, i1 + d1
                if 0 <= j0 < n0 and 0 <= j1 < n1 and active[j1 * n0 + j0]:
                    rr = find(j1 * n0 + j0)
                    if rr != root:
                        out.add(rr)
        return out

    for both_axes in [False, True]:
        stuck = set()
        while True:
            failing = [
                r
                for r in st
                if active[r]
                and not passes(st[r], nmin, bmin)
                and r not in stuck
                and not is_empty(ssum[r], svar[r])
            ]
            if not failing:
                break
            r = min(failing, key=lambda x: (st[x][0] > 0, st[x][1]))
            nb = neighbours(r, both_axes)
            if not nb:
                # e.g. a whole line along the first axis in the first phase; retried in the second phase
                stuck.add(r)
                continue
            trials = [
                (
                    passes(
                        x := group_stats(ssum[r] + ssum[n], svar[r] + svar[n]),
                        nmin,
                        bmin,
                    ),
                    x[1],
                    n,
                )
                for n in nb
            ]
            ok = [t for t in trials if t[0]]
            populated = [m for m in nb if not is_empty(ssum[m], svar[m])]
            if ok:
                n = min(ok, key=lambda t: t[1])[2]
            elif populated:
                n = min(populated, key=lambda m: st[m][1])
            else:
                n = min(nb)  # only empty neighbours: cross them
            a, b = (r, n) if len(members[r]) >= len(members[n]) else (n, r)
            parent[b] = a
            ssum[a] += ssum[b]
            svar[a] += svar[b]
            members[a] += members.pop(b)
            st.pop(b)
            st[a] = group_stats(ssum[a], svar[a])
            stuck.discard(a)
    return np.array([find(c) for c in range(ncell)]).reshape(n1, n0).T


def failing_groups(vals, ulvar, labels, active, nmin, bmin):
    """Number of active, populated groups failing the criteria, recomputed from the bins of every group"""
    ng = labels.max() + 1
    ssum = np.zeros((ng, vals.shape[-1]))
    np.add.at(ssum, labels.ravel(), vals.reshape(-1, vals.shape[-1]))
    svar = np.bincount(labels.ravel(), weights=ulvar.ravel(), minlength=ng)
    gactive = np.bincount(labels.ravel(), weights=active.ravel(), minlength=ng) > 0
    return sum(
        gactive[g]
        and not is_empty(ssum[g], svar[g])
        and not passes(group_stats(ssum[g], svar[g]), nmin, bmin)
        for g in range(ng)
    )


def merge_groups(h, merge_axes, nmin, bmin, ranges=None, helicity_axis="helicity"):
    """Group labels over the non-helicity axes of h, merging only along the two merge_axes and only bins within
    ranges ({axis name: (low, high)}); returns the labels, the axis names, the number of groups that still fail the
    criteria and the mask of the bins within ranges."""
    names = [n for n in h.axes.name if n != helicity_axis]
    hel = h.axes.name.index(helicity_axis)
    vals0 = np.moveaxis(h.values(), hel, -1)
    ulvar0 = np.moveaxis(h.variances(), hel, -1)[..., 0]
    active0 = np.ones(ulvar0.shape, dtype=bool)
    for name, (low, high) in (ranges or {}).items():
        edges = h.axes[name].edges
        shape = [1] * len(names)
        shape[names.index(name)] = -1
        active0 &= ((edges[:-1] >= low) & (edges[1:] <= high)).reshape(shape)
    i0, i1 = names.index(merge_axes[0]), names.index(merge_axes[1])
    others = [i for i in range(len(names)) if i not in (i0, i1)]
    order = [i0, i1, *others]
    vals = np.transpose(vals0, (*order, len(names)))
    ulvar = np.transpose(ulvar0, order)
    active = np.transpose(active0, order)
    labels = np.zeros(ulvar.shape, dtype=np.int64)
    offset = 0
    for idx in itertools.product(
        *[range(ulvar.shape[k]) for k in range(2, ulvar.ndim)]
    ):
        sl = (slice(None), slice(None), *idx)
        lab = merge_plane(vals[sl], ulvar[sl], active[sl], nmin, bmin)
        _, inv = np.unique(lab, return_inverse=True)
        labels[sl] = inv.reshape(lab.shape) + offset
        offset += inv.max() + 1
    labels = np.transpose(labels, np.argsort(order))
    nfail = failing_groups(vals0, ulvar0, labels, active0, nmin, bmin)
    return labels, names, nfail, active0


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--minnloFile",
        required=True,
        help="MiNNLO gen file with the helicity cross sections",
    )
    parser.add_argument(
        "--minnloh", default="nominal_gen", help="histogram with a helicity axis"
    )
    parser.add_argument(
        "--minnloProcs",
        nargs="+",
        default=["Zmumu"],
        help="MiNNLO samples, read as <name>_<era>",
    )
    parser.add_argument("--eras", nargs="+", default=["13TeVGen"])
    parser.add_argument(
        "--noCombineTau", action="store_true", help="Do not add the tau samples"
    )
    parser.add_argument(
        "--mergeAxes",
        nargs=2,
        default=["ptVgen", "massVgen"],
        help="Axes along which bins are merged, the first one first",
    )
    parser.add_argument(
        "--nmin",
        type=float,
        default=100.0,
        help="Minimum effective number of events per group",
    )
    parser.add_argument(
        "--bmin",
        type=float,
        default=0.1,
        help="Minimum of the angular bracket per group",
    )
    parser.add_argument(
        "--massRange",
        nargs=2,
        type=float,
        default=None,
        help="Generated mass range of the samples (massVgen); bins outside are not merged, "
        "e.g. 50 13000 for a reference of the m > 50 GeV samples only",
    )
    parser.add_argument(
        "-o", "--outfile", required=True, help="Output npz file with the group labels"
    )
    parser.add_argument("-v", "--verbose", type=int, default=3, choices=[0, 1, 2, 3, 4])
    args = parser.parse_args()
    logger = logging.setup_logger(__file__, args.verbose)

    h = input_tools.read_mu_hist_combine_tau(
        args.minnloFile,
        args.minnloProcs[0],
        args.minnloh,
        eras=args.eras,
        combine_with_tau=not args.noCombineTau,
    )
    for proc in args.minnloProcs[1:]:
        h = h + input_tools.read_mu_hist_combine_tau(
            args.minnloFile,
            proc,
            args.minnloh,
            eras=args.eras,
            combine_with_tau=not args.noCombineTau,
        )

    ranges = {"massVgen": args.massRange} if args.massRange else None
    labels, names, nfail, active = merge_groups(
        h, args.mergeAxes, args.nmin, args.bmin, ranges
    )
    sizes = np.bincount(labels.ravel())
    merged = sizes[labels] > 1
    logger.info(
        f"{labels.size} bins ({int((~active).sum())} outside the mass range) -> {len(sizes)} groups "
        f"(nmin {args.nmin:g}, bmin {args.bmin:g}), {int(merged.sum())} bins ({100 * merged.mean():.1f}%) in merged "
        f"groups, largest group {sizes.max()} bins, {nfail} groups failing after merging"
    )
    np.savez(
        args.outfile,
        labels=labels,
        axes=np.array(names),
        **{f"edges_{n}": h.axes[n].edges for n in names},
        nmin=args.nmin,
        bmin=args.bmin,
        **({"massRange": args.massRange} if args.massRange else {}),
    )
    logger.info(f"Wrote {args.outfile}")


if __name__ == "__main__":
    main()
