import os
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({
    'font.size': 28,
    # 'font.weight': 'bold',
    'axes.labelsize': 28,
    # 'axes.labelweight': 'bold',
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

# datasets = ["DBLP", "IMDB", "FourSquare"]
datasets = ["DBLP", "IMDB", "FourSquare"]
algorithms = ["L2R", "Naive", "MNC-SAL", "MetaAC-SAL", "DMap-SAL", "LGraph-SAL", "RSE-CSR", "RoseMM"]

# 定义每个算法对应的文件后缀
file_suffix = {
    "L2R": 'l2r_eff',
    "Naive": 'naive_eff',
    "MNC-SAL": 'mnc_eff',
    "MetaAC-SAL": 'meta_eff',
    "DMap-SAL": 'dm_eff',
    "LGraph-SAL": 'lg_eff',
    "RSE-CSR": 'rescsr_eff',
    "RoseMM": 'resmcm_eff'
}

L = [7]

root = os.path.abspath(".")
# print("根目录:", root)

# 假设每个数据集的文件路径遵循相似的结构，例如：
base_fileroot = os.path.join(root, "../../materials/output/efficient")

# 初始化数据字典
data = {algo: [] for algo in algorithms}

for dataset in datasets:
    fileroot = os.path.join(base_fileroot, dataset)
    for algo in algorithms:
        s = 0
        num = 0;
        for l in L:
            file_path = os.path.join(fileroot, file_suffix[algo] + "-l" + str(l) + "-t64.txt")
            timings = readFromTXT(file_path)
            s += sum(timings) * 1000  # 转换为毫秒
            num +=  len(timings)
            num +=  len(timings)
        if(num == 0 ): data[algo].append(None)
        else: data[algo].append(s/num)

# 设置柱状图参数
labels = datasets
num_algorithms = len(algorithms)
x = np.arange(len(labels))  # 数据集数量
bar_width = 0.1  # 每个柱子的宽度
opacity = 0.8

# 颜色列表，可以根据需要调整
colors = ['#001871', '#FF585D', '#FFB549', '#41B6E6', '#629460', '#D0A727', '#5B0888', '#CD6133']
hatches = ['//'    , '-'    , 'xx'     , '++'     , '\\'   , '.'   , '\\'   , '*']

fig, ax = plt.subplots(figsize=(20, 10))

# 为每个算法绘制柱状图
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
plt.ylim((1, 1e7))

ax.set_ylabel('Running time (ms)')

plt.tight_layout(rect=[0, 0, 0.95, 0.95])

# 添加整体标题
fig.text(0.5, 0.01, 'Fig. 12. Efficiency of all SMCM algorithms on all datasets.', ha='center', fontsize=28)

ax.set_xticks(x)
ax.set_xticklabels(labels)

fig.legend(ncol=4, loc='upper center', fontsize=28, frameon=False)

plt.subplots_adjust(left=0.1, right=0.9, top=0.85, bottom=0.15)

root = os.path.abspath(".")
# 保存图像
plt.savefig(root + '/../../materials/output/fig/fig_12_runtime_of_algorithms_on_all_datasets.png', dpi=300)