
import os
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({
    'font.size': 28,
    'axes.labelsize': 28,
    'axes.titlesize': 28,
    'legend.fontsize': 28,
    'xtick.labelsize': 28,
    'ytick.labelsize': 28
})



def readFromTXT(file):
    res = []
    if os.path.exists(file):
        with open(file, 'r') as f:
            for line in f:
                parts = line.split()
                ele = float(parts[2])
                res.append(ele)
    return res

datasets = ["DBLP", "IMDB", "FourSquare"]


algorithms = ["RoseMM", "RSE-CSR", "MKL-CSR", "MKL-CSC", "MKL-BSR"]

file_suffix = {
    "RoseMM": 'resmcm_eff',
    "RSE-CSR": 'rescsr_eff',
    "MKL-CSR": 'mkl_csr_eff',
    "MKL-CSC": 'mkl_csc_eff',
    "MKL-BSR": 'mkl_bsr_eff',
}

L = [5]

root = os.path.abspath(".")

base_fileroot = os.path.join(root, "../../materials/output/efficient")

data = {algo: [] for algo in algorithms}

for dataset in datasets:
    fileroot = os.path.join(base_fileroot, dataset)
    for algo in algorithms:
        s = 0
        num = 0;
        for l in L:
            file_path = os.path.join(fileroot, file_suffix[algo] + "-l" + str(l) + "-t64.txt")
            timings = readFromTXT(file_path)
            s += sum(timings) * 1000 
            num +=  len(timings)
            num +=  len(timings)
        if(num == 0 ): data[algo].append(None)
        else: data[algo].append(s/num)

labels = datasets
num_algorithms = len(algorithms)
x = np.arange(len(labels))  
bar_width = 0.1  
opacity = 0.8

colors = ['#001871', '#FF585D', '#FFB549', '#41B6E6', '#629460', '#D0A727', '#5B0888', '#CD6133']
hatches = ['//', '-', 'xx', '++', '\\\\']

fig, ax = plt.subplots(figsize=(20, 10))

for i, algo in enumerate(algorithms):
    offset = (i - num_algorithms / 2) * bar_width + bar_width / 2
    bars = ax.bar(
        x + offset, 
        data[algo], 
        bar_width, 
        label=algo, 
        color=colors[i % len(colors)], 
        edgecolor='black', 
        alpha=opacity,
        hatch=hatches[i % len(hatches)]
    )

ax.set_yscale('log')
plt.ylim((1, 1e5))

ax.set_ylabel('Running time (ms)')

plt.tight_layout(rect=[0, 0, 0.95, 0.95])

fig.text(0.5, 0.01, 'Fig. 17. Runtime of SMCM with different sparse matrix multiplication on all datasets (length = 5).', ha='center', fontsize=28)

ax.set_xticks(x)
ax.set_xticklabels(labels)

fig.legend(ncol=5, loc='upper center', fontsize=28, frameon=False)

plt.savefig(root + "/../../materials/output/fig/fig_17_runtime_of_smcm_with_differet_spgemm.png")