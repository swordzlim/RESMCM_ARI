import os
import matplotlib.pyplot as plt
import numpy as np
import sys
import textwrap

def readFromTXT(file):
    res=list()
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

mkl_csr = readFromTXT(fileroot + 'mkl_csr_eff.txt')
mkl_csr = sum(mkl_csr)/len(mkl_csr) * 1000

mkl_csc = readFromTXT(fileroot + 'mkl_csc_eff.txt')
mkl_csc = sum(mkl_csc)/len(mkl_csc) * 1000

mkl_bsr = readFromTXT(fileroot + 'mkl_bsr_eff.txt')
mkl_bsr = sum(mkl_bsr)/len(mkl_bsr) * 1000

rescsr = readFromTXT(fileroot + 'rescsr_eff.txt')
rescsr = sum(rescsr)/len(rescsr) * 1000

resmcm = readFromTXT(fileroot + 'resmcm_eff.txt')
resmcm = sum(resmcm)/len(resmcm) * 1000

data = {
    'RoseMM': resmcm,
    'RSE-CSR': rescsr,
    'MKL-CSR': mkl_csr,
    'MKL-CSC': mkl_csc,
    'MKL-BSR': mkl_bsr,
}

labels = list(data.keys())
values = list(data.values())
x = np.arange(len(labels))
bar_width = 0.35

fig, ax = plt.subplots(figsize=(10, 7))

ax.bar(x, values, bar_width, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'], edgecolor='black')

ax.set_yscale('log')

ax.set_ylabel('Running time (ms)', fontsize=14)
ax.set_xlabel('Algorithms', fontsize=14)
len_arg=sys.argv[1]
title='Runtime of SMCM with different sparse matrix multiplication on all datasets with length=' + len_arg
wrapped_title = "\n".join(textwrap.wrap(title, width=60))
ax.set_title(wrapped_title, fontsize=16)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=12)


plt.tight_layout()
plt.show()

plt.savefig(root + "/../../materials/output/fig/IMDB/mkl.png")