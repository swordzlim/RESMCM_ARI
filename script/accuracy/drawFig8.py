import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

def readFromTXT(file):
    """Reads the second column of a TXT file and returns a list of floats."""
    res = []
    if os.path.exists(file):
        with open(file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        ele = float(parts[1])
                        res.append(ele)
                    except ValueError:
                        continue  
    else:
        print(f"File not found: {file}")
    return res

plt.rcParams.update({
    'font.size': 28,
    'axes.labelsize': 28,
    'axes.titlesize': 28,
    'legend.fontsize': 28,
    'xtick.labelsize': 28,
    'ytick.labelsize': 28
})

datasets = ['DBLP', 'IMDB', 'FourSquare']
estimators = ['RS-estimator', 'MetaAC', 'MNC', 'DMap', 'LGraph']
colors = ['#629460', '#FF595E', '#41B6E6', '#FFFFA5', '#FFA5FF']

root = "/home/chunxu/cxx/RoseMM_v2/RESMCM_ARI/script/accuracy"

data = {dataset: {estimator: [] for estimator in estimators} for dataset in datasets}

for dataset in datasets:
    accuracy_fileroot_template = os.path.join(root, "../../materials/output/accuracy/{dataset}/")
    efficient_fileroot_template = os.path.join(root, "../../materials/output/efficient/{dataset}/")
    acc_fileroot = accuracy_fileroot_template.format(dataset=dataset)
    gb_fileroot = efficient_fileroot_template.format(dataset=dataset)

    length = 2
    gb_file = os.path.join(gb_fileroot, f'resmcm_eff-l{length}-t64.txt')
    print(f"Ground truth file is: {gb_file}")
    acc_file = os.path.join(acc_fileroot, f'res_acc-l{length}.txt')
    meta_file = os.path.join(acc_fileroot, f'meta_acc-l{length}.txt')
    mnc_file = os.path.join(acc_fileroot, f'mnc_acc-l{length}.txt')
    dm_file = os.path.join(acc_fileroot, f'dm_acc-l{length}.txt')
    lg_file = os.path.join(acc_fileroot, f'lg_acc-l{length}.txt')

    gb = readFromTXT(gb_file)
    if not gb:
        print(f"No baseline data for {dataset} at length {length}. Skipping.")
        continue

    res = readFromTXT(acc_file)
    meta = readFromTXT(meta_file)
    mnc = readFromTXT(mnc_file)
    dm = readFromTXT(dm_file)
    lg = readFromTXT(lg_file)

    min_len = min(len(gb), len(res), len(meta), len(mnc), len(dm), len(lg))
    if min_len == 0:
        print(f"Insufficient data for {dataset} at length {length}. Skipping.")
        continue

    gb = gb[:min_len]
    res = res[:min_len]
    meta = meta[:min_len]
    mnc = mnc[:min_len]
    dm = dm[:min_len]
    lg = lg[:min_len]

    def compute_rel_error(estimator):
        lst = []
        for i in range(len(gb)):
            if estimator[i] == 0:
                lst.append(1e7)
            else:
                val = max(gb[i], estimator[i]) / min(gb[i], estimator[i])
                if val > 1e7:
                    val = 1e7
                lst.append(val)
        return lst

    data[dataset]['RS-estimator'] = sorted(compute_rel_error(res))
    data[dataset]['MetaAC'] = sorted(compute_rel_error(meta))
    data[dataset]['MNC'] = sorted(compute_rel_error(mnc))
    data[dataset]['DMap'] = sorted(compute_rel_error(dm))
    data[dataset]['LGraph'] = sorted(compute_rel_error(lg))

boxplot_data = []
boxplot_positions = []
width = 0.15
offsets = np.linspace(-width*2, width*2, len(estimators))  

for dataset_idx, dataset in enumerate(datasets):
    for est_idx, estimator in enumerate(estimators):
        if len(data[dataset][estimator]) > 0:
            boxplot_data.append(data[dataset][estimator])
            boxplot_positions.append(dataset_idx + offsets[est_idx])

fig, ax = plt.subplots(figsize=(20, 10))
bp = ax.boxplot(boxplot_data, positions=boxplot_positions, widths=width, whis=(0, 100), patch_artist=True, showfliers=True, medianprops=dict(color='black'))

for i, patch in enumerate(bp['boxes']):
    color = colors[i % len(estimators)]
    patch.set_facecolor(color)

ax.set_yscale('log')
ax.set_ylim(0.5, 1e7)

ax.set_ylabel('Relative Error')
ax.set_xticks(range(len(datasets)))
ax.set_xticklabels(datasets, fontsize=28)

ax.grid(True, ls="--", linewidth=0.5)

legend_handles = [mpatches.Patch(color=colors[i], label=estimators[i]) for i in range(len(estimators))]
fig.legend(handles=legend_handles, labels=estimators, ncol=5, fontsize=28, frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1))

plt.tight_layout()

fig.text(0.5, 0.02, 'Fig. 8. Relative error of all estimators for sparse matrix-matrix multiplication on all datasets.', ha='center', fontsize=28)

plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15)

root = os.path.abspath(".")
plt.savefig(root + '/../../materials/output/fig/fig_8_accuracy_comparison_for_smm.png', dpi=300)