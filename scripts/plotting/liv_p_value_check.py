import re

import numpy as np
from scipy import stats

from rabbit import io_tools

# this_file = "/ceph/submit/data/user/j/jbenke/fitresults/fitresults_cd_rapidity_randTime.hdf5"
this_file = "fitresults_dtu_rapidity_randTime.hdf5"

fitresult = io_tools.get_fitresult(this_file)
parms = fitresult["parms"].get()
cov = fitresult["cov"].get()

# Wilson coefficients are named like cxxu, dxys, ctxd, ...; everything else is a nuisance
names = list(parms.axes[0])
wilson = [n for n in names if re.fullmatch(r"[cd]t?(xx|xy|xz|yz|x|y)[uds]", n)]
idx = [names.index(n) for n in wilson]

p = parms.values()[idx]
C = cov.values()[np.ix_(idx, idx)]

chi2 = p @ np.linalg.solve(C, p)
pval = stats.chi2.sf(chi2, df=len(p))

print(f"coefficients: {wilson}")
print(f"values:       {p}")
print(f"uncertainties:{np.sqrt(np.diag(C))}")
print(f"chi2 = {chi2:.3f}, ndf = {len(p)}, p-value = {pval:.4g}")
