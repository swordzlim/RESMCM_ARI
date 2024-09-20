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
                ele = float(parts[2])
                res.append(ele)
    return res


# ("dm" "lg" "meta" "mnc" "rescsr" "resmcm")
root = os.path.abspath(".")
print(root)
fileroot=root + "/../../materials/output/efficient/IMDB/"

l2r = readFromTXT(fileroot + 'l2r_eff.txt')
l2r = sum(l2r)/len(l2r) * 1000

naive = readFromTXT(fileroot + 'naive_eff.txt')
naive = sum(naive)/len(naive) * 1000

mnc = readFromTXT(fileroot + 'mnc_eff.txt')
mnc = sum(mnc)/len(mnc) * 1000

meta = readFromTXT(fileroot + 'meta_eff.txt')
meta = sum(meta)/len(meta) * 1000

dm = readFromTXT(fileroot + 'dm_eff.txt')
dm = sum(dm)/len(dm) * 1000

lg = readFromTXT(fileroot + 'lg_eff.txt')
lg = sum(lg)/len(lg) * 1000

rescsr = readFromTXT(fileroot + 'rescsr_eff.txt')
rescsr = sum(rescsr)/len(rescsr) * 1000

resmcm = readFromTXT(fileroot + 'resmcm_eff.txt')
resmcm = sum(resmcm)/len(resmcm) * 1000

data = {
    'L2R': l2r,
    'Naive': naive,
    'MNC-SAL': mnc,
    'MetaAC-SAL': meta,
    'DMap-SAL': dm,
    'LGraph-SAL': lg,
    'RSE-CSR': rescsr,
    'RoseMM': resmcm,
}

labels = list(data.keys())
values = list(data.values())
x = np.arange(len(labels))
bar_width = 0.35

fig, ax = plt.subplots(figsize=(10, 7))

ax.bar(x, values, bar_width, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'], edgecolor='black')

ax.set_yscale('log')
plt.ylim((1, 10e7))

ax.set_ylabel('Running time (ms)', fontsize=14)
ax.set_xlabel('Algorithms', fontsize=14)
len_arg=sys.argv[1]

title='Runtime of all estimators for SMCM on IMDB dataset with length=' + len_arg
wrapped_title = "\n".join(textwrap.wrap(title, width=60))
ax.set_title(wrapped_title, fontsize=16)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=12)


plt.tight_layout()
plt.show()

plt.savefig(root + "/../../materials/output/fig/IMDB/eff.png")