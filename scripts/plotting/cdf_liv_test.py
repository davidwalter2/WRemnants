import h5py
import hist
import numpy as np

from wremnants.utilities.io_tools import input_tools


def get_cdf(hist_in):
    arr = hist_in.copy().to_numpy()
    data_arr = arr[0]
    mass_data = np.sum(data_arr, axis=0)
    cdf_arr = np.cumsum(mass_data)
    cdf_arr /= cdf_arr[-1]
    return cdf_arr, data_arr.shape


indir_data = "/work/submit/jbenke/WRemnants/scripts/histmakers/"
infile_data = indir_data + "mz_dilepton_liv_scetlib_dyturboCorr.hdf5"

h5file = h5py.File(infile_data, "r")
results = input_tools.load_results_h5py(h5file)
hist_in = results["ZmumuPostVFP"]["output"]["fine_bin_axis_gen"].get()


def get_cdf(hist_in):
    arr = hist_in.copy().to_numpy()
    data_arr = arr[0]
    cdf_arr = np.cumsum(data_arr)
    cdf_arr /= cdf_arr[-1]

    return cdf_arr, data_arr.shape


def make_quantiles(hist_in, n_quantiles, axis_name):
    cdf_output, hist_shape = get_cdf(hist_in)

    if ((hist_shape[0]) % n_quantiles) != 0:
        print(
            "wrong number of quantiles. choose something that factors into %s"
            % hist_shape[0]
        )
    else:
        cdf_vals_in = np.linspace(0, 1, n_quantiles + 1)
        x_vals = hist_in.axes[0].edges[1:]
        new_edges = np.interp(cdf_vals_in, cdf_output, x_vals)
        new_edges = [*new_edges]
        new_axis = hist.axis.Variable(new_edges, name=axis_name)
        return new_axis


# pdb.set_trace()
k = make_quantiles(hist_in, 6, "pt_sublead")
print(k)


## weighted average for pt on 9/12/25
a = [
    25.1375,
    27.3689,
    29.5599,
    31.7123,
    33.8092,
    35.8465,
    37.8381,
    39.7941,
    41.7513,
    43.8944,
    80,
]  # dtst
b = [
    24.9125,
    34.9961,
    38.9813,
    41.8262,
    43.9753,
    45.6349,
    47.1672,
    49.1107,
    52.1681,
    58.071,
    80,
]  # dtdt
c = [
    25.1375,
    27.2737,
    29.4016,
    31.5117,
    33.5912,
    35.6305,
    37.6446,
    39.6377,
    41.6384,
    43.8174,
    80,
]  ## stst
avg = []
b_weight = 8753341.531365383
a_weight = 2071457.4297493682
c_weight = 2519865.2781415316
weightsum = a_weight + b_weight + c_weight


a_weight /= weightsum
b_weight /= weightsum
c_weight /= weightsum
for i in range(len(a)):
    avg.append(round((a[i] * a_weight + b[i] * b_weight + c[i] * c_weight), 5))
# print(avg)
### new pt bins:  [15, 24, 32.35393, 35.70991, 38.30856, 40.43642, 42.22635, 43.92092, 45.87573, 48.56281, 53.1789, 80]

a = [-2.4, -1.63869, -0.849581, -0.0505165, 0.773543, 1.58705, 2.4]  ## dtst
c = [-2.4, -1.57806, -0.802603, -0.0277886, 0.740474, 1.53275, 2.4]  ## stst

b = [-2.4, -1.30224, -0.606947, 0.00701935, 0.622103, 1.31844, 2.4]  ##dtdt
