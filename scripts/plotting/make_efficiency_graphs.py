import pickle
import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.rcParams.update({"font.size": 12})


with open("efficiency_values.pkl", "rb") as f:
    input_dict = pickle.load(f)
### so these get loaded in as pt, eta.

date = "2025-10-08"
postfit = False  ### I SHOULD DO THIS
print(date)

# pt_cutoff = 2
pt_cutoff = 0

# Read your input from a text file
with open("efficiency_graph_fit_9-5-25.txt", "r") as mult:
    input_text = mult.read()

# Extract parameter names and constraints using regular expressions
pattern = re.compile(
    r"^\s*(\S+)\s+([-+]?\d*\.?\d+)\s+\+/-\s+([-+]?\d*\.?\d+)", re.MULTILINE
)
matches = pattern.findall(input_text)
parameters = [param for param, _, _ in matches]
param_constraints = [constraint for _, _, constraint in matches]

params_trig = [param for param in parameters if "trig_" in param]
params_n = [param for param in parameters if "n_" in param]
params_id = [param for param in parameters if "id_" in param]


def extract_ijk(text):
    # Match digits after 'mll' and 'time'
    pt_match = re.search(r"pt(\d+)", text)
    eta_match = re.search(r"eta(\d+)", text)

    time_match = re.search(r"time(\d+)", text)

    pt = int(pt_match.group(1)) if pt_match else None
    eta = int(eta_match.group(1)) if eta_match else None

    time = int(time_match.group(1)) if time_match else None
    return pt, eta, time


# eps_prime = np.zeros(
#     [
#         len(input_dict.keys()),
#         input_dict["stst_id"].shape[0],
#         input_dict["stst_id"].shape[1],
#         input_dict["stst_id"].shape[2],
#     ]
# )
# all_param_list = [params_n, params_trig, params_id]

# for k in range(len(all_param_list)):
#     for name in all_param_list[k]:
#         i = parameters.index(name)
#         value = param_constraints[i]
#         pt, eta, time = extract_ijk(name)
#         eps_prime[k, pt, eta, time] = float(value)
#         eps_prime[k + 3, pt, eta, time] = float(value)
#         eps_prime[k + 6, pt, eta, time] = float(value)

##### looks at a sing
key_list = list(input_dict.keys())


##### this is all prefit


pt_bins = [
    15,
    25,
    32.35393,
    35.70991,
    38.30856,
    40.43642,
    42.22635,
    43.92092,
    45.87573,
    48.56281,
    53.1789,
    80,
]  # , 80,], ### these are not the pt bins

pt_bins = [
    15,
    # 21,
    25,
    45,
    65,
    80,
]

num_eta_bins = 3
# pdb.set_trace()
pt_bins = [(pt_bins[i] + (pt_bins[i + 1])) / 2 for i in range(len(pt_bins) - 1)]
eta_bins = np.linspace(-2.4, 2.4, num_eta_bins)

xx, yy = np.meshgrid(eta_bins, pt_bins)
zz = np.zeros((len(pt_bins), num_eta_bins))
print(zz.shape)
for efficiency in range(len(key_list)):
    for key in [key_list[efficiency]]:
        # pdb.set_trace()
        # if "epsilon" not in key:
        #     for j in range(len(input_dict[key][0][0])):  # eta bin
        #         temp = []
        #         for i in range(len(input_dict[key][0])):  ## pt bin
        #             temp.append(np.max(input_dict[key][:, i, j]))
        #         if "high" in key:
        #             plt.plot(
        #                 pt_bins[pt_cutoff:], temp[pt_cutoff:], label=f"eta bin: {j}", color=f"C{j}"
        #             )
        #         else:
        #             plt.plot(pt_bins, temp, label=f"eta bin: {j}", color=f"C{j}")
        #     plt.title(key + " prefit")
        #     plt.legend()
        #     # plt.ylim([0.99, 1.01])
        #     plt.xlabel("pt bin lower edge [GeV]")
        #     plt.ylabel("input efficiency")
        #     plt.tight_layout()
        #     plt.savefig(
        #         f"/home/submit/jbenke/public_html/liv_uncert/efficiency/{date}/pt_projection/eff_{key}_pt_projection.png"
        #     )
        #     plt.clf()

        if "epsilon" in key:
            for j in range(len(input_dict[key][0][0])):  # eta bin
                temp = []
                for i in range(len(input_dict[key][0])):  ## pt bin
                    temp.append(np.average(input_dict[key][:, i, j]))
                # pdb.set_trace()
                if (
                    "high" in key
                ):  ###SHOULD COME UP WITH A SMART WAY TO PLOT THESE SIMULTANEOUSLY
                    plt.plot(
                        pt_bins[pt_cutoff:],
                        temp[pt_cutoff:],
                        label=f"eta bin: {j}",
                        color=f"C{j}",
                    )
                if "low" in key:
                    plt.plot(
                        pt_bins[:pt_cutoff],
                        temp[:pt_cutoff],
                        label=f"eta bin: {j}",
                        color=f"C{j}",
                    )
            plt.title(key + " prefit")
            plt.legend()
            plt.xlabel("pt bin lower edge [GeV]")
            plt.ylabel("input efficiency")
            plt.tight_layout()
            plt.savefig(
                f"/home/submit/jbenke/public_html/liv_uncert/efficiency/{date}/pt_projection/{key}.png"
            )
            plt.clf()

        if "true" in key or "COMBINED" in key:
            for j in range(len(input_dict[key][0][0])):  # eta bin
                temp = []
                for i in range(len(input_dict[key][0])):  ## pt bin
                    temp.append(np.average(input_dict[key][:, i, j]))
                    zz[i, j] = np.average(input_dict[key][:, i, j])

                plt.plot(pt_bins[:], temp[:], label=f"eta bin: {j}", color=f"C{j}")

            plt.title(key + " prefit")
            plt.legend()
            # plt.ylim([0.99, 1.01])
            plt.xlabel("pt bin lower edge [GeV]")
            plt.ylabel("input efficiency")
            plt.tight_layout()
            plt.savefig(
                f"/home/submit/jbenke/public_html/liv_uncert/efficiency/{date}/pt_projection/{key}.png"
            )
            plt.clf()


# pdb.set_trace()
plt.clf

# plt.plot(xx, yy, marker = 'o')
# heatmap = np.array([xx, yy, zz])
df = pd.DataFrame(zz, columns=eta_bins, index=pt_bins)
ax = sns.heatmap(df)
ax.invert_yaxis()
# plt.contourf(eta_bins, pt_bins, zz)
plt.xlabel("eta")
plt.ylabel("pt")
# plt.colorbar()
plt.tight_layout()
plt.savefig(
    f"/home/submit/jbenke/public_html/liv_uncert/efficiency/{date}/pt_projection/efficiency_heatmap.png"
)


for efficiency in range(len(key_list)):
    for key in [key_list[efficiency]]:
        if "epsilon" in key:
            for i in range(len(input_dict[key][0])):  ## pt bin
                temp = []
                for j in range(len(input_dict[key][0][0])):  # eta bin
                    temp.append(np.average(input_dict[key][:, i, j]))
                plt.plot(eta_bins[:], temp[:], label=f"pt bin: {i}", color=f"C{i}")
            plt.title(key + " prefit")
            plt.legend()
            # plt.ylim([0.99, 1.01])
            plt.xlabel("eta bin lower edge")
            plt.ylabel("input efficiency")
            plt.tight_layout()
            plt.savefig(
                f"/home/submit/jbenke/public_html/liv_uncert/efficiency/{date}/eta_projection/{key}.png"
            )
            plt.clf()
        if "true" in key:
            for i in range(len(input_dict[key][0])):  ## pt bin
                temp = []
                for j in range(len(input_dict[key][0][0])):  # eta bin
                    temp.append(np.max(input_dict[key][:, i, j]))
                plt.plot(eta_bins, temp, label=f"pt bin: {i}", color=f"C{i}")

            plt.title(key + " prefit")
            plt.legend()
            # plt.ylim([0.99, 1.01])
            plt.xlabel("eta bin lower edge")
            plt.ylabel("input efficiency")
            plt.tight_layout()
            plt.savefig(
                f"/home/submit/jbenke/public_html/liv_uncert/efficiency/{date}/eta_projection/{key}.png"
            )
            plt.clf()


"""
if postfit:
    #########################################333
    ##### this is all postfit

    for efficiency in range(len(key_list)):

        for key in [key_list[efficiency]]:
            for j in range(0, input_dict[key].shape[1]):  # eta bin
                temp = []
                for i in range(2, input_dict[key].shape[0]):  ## pt bin
                    temp.append(
                        np.max(
                            (input_dict[key][i, j] - 1) * eps_prime[efficiency, i, j, :]
                            + 1
                        )
                    )
                    # pdb.set_trace()
                plt.plot(pt_bins, temp, label=f"eta bin: {j}", color=f"C{j}")

        plt.title(key + " prefit")
        plt.legend()
        # plt.ylim([0.99, 1.01])
        plt.xlabel("pt bin lower edge [GeV]")
        plt.ylabel("input efficiency")
        plt.tight_layout()
        plt.savefig(
            f"/home/submit/jbenke/public_html/liv_uncert/efficiency/{date}/pt_projection/postfit/eff_{key}_pt_projection.png"
        )
        plt.clf()

    ####
"""

"""
# Read your input from a text file
with open("efficiency_graph_fit_9-4-25.txt", "r") as f:
    input_text = f.read()

# Extract parameter names and constraints using regular expressions
pattern = re.compile(
    r"^\s*(\S+)\s+([-+]?\d*\.?\d+)\s+\+/-\s+([-+]?\d*\.?\d+)", re.MULTILINE
)
matches = pattern.findall(input_text)
parameters = [param for param, _, _ in matches]
param_constraints = [constraint for _, _, constraint in matches]

params_hlt = [param for param in parameters if "hlt_" in param]
params_id = [param for param in parameters if "id_" in param]
params_num = [param for param in parameters if "n_" in param]


def extract_mll_and_time(text):
    # Match digits after 'mll' and 'time'
    pt_match = re.search(r"pt(\d+)", text)
    eta_match = re.search(r"eta(\d+)", text)

    time_match = re.search(r"time(\d+)", text)

    pt = int(pt_match.group(1)) if pt_match else None
    eta = int(pt_match.group(1)) if eta_match else None

    time = int(time_match.group(1)) if time_match else None
    return mll, time


colors = ["C0", "C1", "C2"]
check_color = [0, 0, 0]
legend_entries = ["m = 89.3 - 90.1 GeV", "m = 90.1 - 90.8 GeV", "m = 90.8 - 91.4 GeV"]
low_mll_hlt = []
med_mll_hlt = []
high_mll_hlt = []
for name in params_hlt:
    i = parameters.index(name)
    value = param_constraints[i]
    mll, time = extract_mll_and_time(name)

    if mll == 3:
        low_mll_hlt.append(float(value) * 0.01)
    elif mll == 4:
        med_mll_hlt.append(float(value) * 0.01)
    elif mll == 5:
        high_mll_hlt.append(float(value) * 0.01)

plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

plt.plot(low_mll_hlt, color="C0", label="m = 89.3 - 90.1 GeV")
plt.plot(med_mll_hlt, color="C1", label="m = 90.1 - 90.8 GeV")
plt.plot(high_mll_hlt, color="C2", label="m = 90.8 - 91.4 GeV")

plt.xlabel("Sidereal time [hr]")
plt.ylabel("HLT efficiency uncertainty")
plt.legend(loc="center right")
plt.tight_layout()

plt.savefig("/home/submit/jbenke/public_html/liv_uncert/efficiency/hlt_efficencies.png")
plt.clf()


low_mll_id = []
med_mll_id = []
high_mll_id = []
for name in params_id:
    i = parameters.index(name)
    value = param_constraints[i]
    mll, time = extract_mll_and_time(name)

    if mll == 3:
        low_mll_id.append(float(value) * 0.01)
    elif mll == 4:
        med_mll_id.append(float(value) * 0.01)
    elif mll == 5:
        high_mll_id.append(float(value) * 0.01)

plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

plt.plot(low_mll_id, color="C0", label="m = 89.3 - 90.1 GeV")
plt.plot(med_mll_id, color="C1", label="m = 90.1 - 90.8 GeV")
plt.plot(high_mll_id, color="C2", label="m = 90.8 - 91.4 GeV")

plt.xlabel("Sidereal time [hr]")
plt.ylabel("ID efficiency uncertainty")
plt.legend(loc="center right")
plt.tight_layout()

plt.savefig("/home/submit/jbenke/public_html/liv_uncert/efficiency/id_efficencies.png")


# print(output)

# Prepare CSV data
# csv_data = [("Parameter", "Constraint")]
# csv_data.extend([(param, constraint) for param, _, constraint in matches])

# # Save to a CSV file
# with open("efficiency_output.csv", "w", newline="") as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerows(csv_data)

"""
