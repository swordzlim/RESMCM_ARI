import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    'font.size': 28,
    'axes.labelsize': 28,
    'axes.titlesize': 28,
    'legend.fontsize': 28,
    'xtick.labelsize': 28,
    'ytick.labelsize': 28
})

colors = {
    'L2R': '#001871',          # c1
    'Naive': '#FF585D',        # c2
    'MNC-SAL': '#FFB549',      # c3
    'MetaAC-SAL': '#41B6E6',   # c4
    'DMap-SAL': '#629460',     # c5
    'LGraph-SAL': '#D0A727',   # c6
    'RSE-CSR': '#5B0888',      # c7
    'RoseMM': '#CD6133'         # c8
}

markers = {
    'L2R': 'x',
    'Naive': '^',
    'MNC-SAL': 'D',
    'MetaAC-SAL': 's',
    'DMap-SAL': 'P',
    'LGraph-SAL': 'v',
    'RSE-CSR': 'o',
    'RoseMM': 'p'
}

def readFromTXT(file):
    lst=list()
    if os.path.exists(file):
        with open(file, 'r') as f:
            for line in f:
                parts = line.split()
                ele = float(parts[2])
                lst.append(ele)
    return lst


def get_data(ds):
    root = os.path.abspath(".")
    eff_fileroot=root + "/../../materials/output/efficient/" + ds +"/"
    x = [8, 16, 32, 64]

    algs = ['l2r', 'naive', 'mnc', 'meta', 'dm', 'lg', 'rescsr', 'resmcm']
    array2d = list()
    for alg in algs:
        array2d.append(list())
        for t in x:
            iFile = eff_fileroot + alg + '_eff-l7-t' + str(t) + '.txt'
            print(iFile)
            lst = readFromTXT(iFile)
            sum = 0
            for item in lst:
                sum += item
            if(len(lst) == 0):
                mean = None
            else:
                mean = sum / len(lst)
                mean = mean * 1000
            array2d[-1].append(mean)

    data = {
        'size': x,
        'l2r': array2d[0],
        'naive': array2d[1],
        'MNC2D': array2d[2],
        'META2D': array2d[3],
        'DMap2D': array2d[4],
        'LGraph2D': array2d[5],
        'RES-CSR': array2d[6],
        'RES2D': array2d[7]
    }
    return data

dataset_names = ['DBLP', 'IMDB', 'FourSquare']

datasets = {
    ds: get_data(ds) for ds in dataset_names
}

fig, axes = plt.subplots(1, 3, figsize=(40, 10))

num = 0
for ax, name in zip(axes, dataset_names):
    data = datasets[name]
    x = data['size']
    
    for method in colors.keys():
        if method == 'RoseMM':
            label = 'RoseMM' 
            y = data['RES2D']
        elif method == 'LGraph-SAL':
            label = 'LGraph-SAL'
            y = data['LGraph2D']
        elif method == 'DMap-SAL':
            label = 'DMap-SAL'
            y = data['DMap2D']
        elif method == 'MNC-SAL':
            label = 'MNC-SAL'
            y = data['MNC2D']
        elif method == 'MetaAC-SAL':
            label = 'MetaAC-SAL'
            y = data['META2D']
        elif method == 'RSE-CSR':
            label = 'RSE-CSR'
            y = data['RES-CSR']
        elif method == 'Naive':
            label = 'Naive'
            y = data['naive']
        elif method == 'L2R':
            label = 'L2R'
            y = data['l2r']
        
        ax.plot(x, y, marker=markers[method],
                color=colors[method], linewidth=4.5, markersize=20, label=label)
    
    ax.set_yscale('log')
    ax.set_xlim(6, 66)
    ax.set_ylabel('Running time (ms)')
    ax.set_title("(" + chr(ord("a") + num) + ") "+ name, y=-0.24)
    num += 1
    ax.grid(True, which="both", ls="--", linewidth=0.5)

    ax.set_xlabel('Length')

for ax in axes:
    ax.set_xticks([8, 16, 32, 64])
    ax.set_xticklabels([8, 16, 32, 64], fontsize=28)

handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=4, fontsize=28, frameon=False)

plt.tight_layout(rect=[0, 0, 0.95, 0.95])

fig.subplots_adjust(hspace=0.1, top = 0.85, bottom=0.25)

fig.text(0.5, 0.05, 'Fig. 14. Effect of the number of threads.', ha='center')


root = os.path.abspath(".")
plt.savefig(root + '/../../materials/output/fig/fig_14_matrix_chain_thread_effect.png', dpi=300)