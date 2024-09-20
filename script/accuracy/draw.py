import os
import matplotlib.pyplot as plt
import numpy as np
import sys
import textwrap

def readFromTXT(file):
    res=list()
    if os.path.exists(file):
        with open(file, 'r') as f:
            for line in f:
                parts = line.split()
                ele = float(parts[1])
                res.append(ele)
    return res


# ("dm" "lg" "meta" "mnc" "rescsr" "resmcm")
root = os.path.abspath(".")
print(root)
acc_fileroot=root + "/../../materials/output/accuracy/IMDB/"
gb_fileroot=root + "/../../materials/output/efficient/IMDB/"

gb=readFromTXT(gb_fileroot + 'resmcm_eff.txt')

res = readFromTXT(acc_fileroot + 'res_acc.txt')
res = [max(gb[i],res[i])/min(gb[i],res[i]) for i in range(len(res))]

meta = readFromTXT(acc_fileroot + 'meta_acc.txt')
meta = [max(gb[i],meta[i])/min(gb[i],meta[i]) for i in range(len(meta))]

mnc = readFromTXT(acc_fileroot + 'mnc_acc.txt')
mnc = [max(gb[i],mnc[i])/min(gb[i],mnc[i]) for i in range(len(mnc))]

dm = readFromTXT(acc_fileroot + 'dm_acc.txt')
dm = [max(gb[i],dm[i])/min(gb[i],dm[i]) for i in range(len(dm))]

lg = readFromTXT(acc_fileroot + 'lg_acc.txt')
lg = [max(gb[i],lg[i])/min(gb[i],lg[i]) for i in range(len(lg))]




data = {
    'RS-estimator': res,
    'MetaAC': meta,
    'MNC': mnc,
    'DMap': dm,
    'LGraph': lg,
}

labels = list(data.keys())
values = list(data.values())
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

fig, ax = plt.subplots(figsize=(10, 7))

boxplots = ax.boxplot(values, patch_artist=True, showfliers=False, 
                      medianprops=dict(color='red'))

for patch, color in zip(boxplots['boxes'], colors):
    patch.set_facecolor(color)

ax.set_yscale('log')
plt.ylim((1, 10e7))

ax.set_ylabel('Running time (ms)', fontsize=14)
ax.set_xlabel('Algorithms', fontsize=14)

len_arg = sys.argv[1]

title='Relative error of all estimators for SMCM on IMDB dataset with length=' + len_arg
wrapped_title = "\n".join(textwrap.wrap(title, width=60))
ax.set_title(wrapped_title, fontsize=16)

ax.set_xticks(np.arange(1, len(labels) + 1))
ax.set_xticklabels(labels, fontsize=12)

plt.tight_layout()
plt.show()

plt.savefig(root + "/../../materials/output/fig/IMDB/acc.png")