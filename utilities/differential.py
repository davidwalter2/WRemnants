from utilities import logging
import hist

logger = logging.child_logger(__name__)

axis_xnorm = hist.axis.Regular(1, 0., 1., name = "count", underflow=False, overflow=False)

def get_pt_eta_axes(n_bins, min_pt, max_pt, max_eta):

    # gen axes for differential measurement
    if len(n_bins) == 2:
        genBinsPt = n_bins[0]
        genBinsEta = n_bins[1]
    else:
        raise IOError(f"Specified format 'n_bins {n_bins}' can not be processed! Please specify the number of gen bins for pT and |eta|")
    axis_ptGen = hist.axis.Regular(genBinsPt, min_pt, max_pt, name = "ptGen")
    axis_etaGen = hist.axis.Regular(genBinsEta, 0, max_eta, underflow=False, name = "etaGen")
    axes = [axis_ptGen, axis_etaGen]
    cols = ["ptGen", "etaGen"]

    return axes, cols

def get_pt_eta_charge_axes(n_bins, min_pt, max_pt, max_eta):

    axes, cols = get_pt_eta_axes(n_bins, min_pt, max_pt, max_eta)

    axis_qGen = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "qGen")

    axes.append(axis_qGen)
    cols.append("qGen")

    return axes, cols

# other gen axes
def get_gen_axes(n_bins, gen_vars):

    # gen axes for differential measurement
    if isinstance(n_bins, int) and isinstance(gen_vars,str):
        n_bins = [n_bins]
        gen_vars = [gen_vars]
    elif len(n_bins) == len(gen_vars):
        pass
    else:
        raise IOError(f"Specified format 'n_bins {n_bins}' and 'axes' {gen_vars} can not be processed! Please specify the numbers of gen bins for given gen axes")

    cols = gen_vars
    axes = [make_gen_axis(v,n) for v,n in zip(gen_vars, n_bins)]

    return axes, cols

def make_gen_axis(name, n_bins=None, bin_lo=None, bin_hi=None):
    gen_axes = {
        "mVGen": {
            "n_bins": 60,
            "bin_lo": 60,
            "bin_hi": 120
        },
        "absYVGen": {
            "n_bins": 25,
            "bin_lo": 0,
            "bin_hi": 2.5,
            "underflow": False
        }, 
        "ptVGen": {
            "n_bins": 44,
            "bin_lo": 26,
            "bin_hi": 70
        }    
    }

    info = gen_axes[name]
    n_bins = n_bins if n_bins else info["bins"]
    bin_lo = bin_lo if bin_lo else info["bin_lo"]
    bin_hi = bin_hi if bin_hi else info["bin_hi"]
    underflow = info.get("underflow", True)
    overflow = info.get("overflow", True)

    return hist.axis.Regular(n_bins, bin_lo, bin_hi, underflow=underflow, overflow=overflow, name=name)
