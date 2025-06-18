import os
import numpy as np
import pandas as pd
from matplotlib.font_manager import fontManager
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mpt
import matplotlib.cm as mpm
from itertools import product
import seaborn as sns
import sys
import argparse
import warnings

warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", category=UserWarning)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run ATAC sample replicate correlation")
    parser.add_argument("--homer_results",required=True, type=str, help="The homer results directory")
    parser.add_argument("--group",required=True,type=str, help="The grouping file of homer results")
    parser.add_argument("--motifs", required=True,type=str, help="The text file to specify the motifs to plot")
    parser.add_argument("-o",required=True,type=str, help="The output directory")
    parser.add_argument("--log_p_filter",default=-2,type=float, help="The filtering log p value")
    parser.add_argument("--log_p_point",default=-5,type=float, help="The log p threshold to plot stars")
    parser.add_argument("--log_p_min",default=-50,type=float, help="The minimum log p value to plot")
    parser.add_argument("--ratio",default=3,type=float, help="The ratio threshold to plot stars")
    return parser.parse_args()

args = parse_arguments()
motif_dir = args.homer_results
save_dir = args.o
group_file = args.group
choose_motifs = args.motifs #required,，所有的TF会产出太长的图而报错
log_p_value = args.log_p_filter #过滤
log_p_value_point = args.log_p_point # 画星星
log_p_value_min= args.log_p_min #画图
ratio_min= args.ratio # 画星星

print("Plot motif: {}".format(motif_dir))
data = list()

group = pd.read_csv(group_file, sep=None, engine='python',)

name_handle = lambda d: d.replace(")", "").replace("?", "").replace(",,", ",").strip(",").strip()

missing_elements = set(group['dir']) - set(os.listdir(motif_dir))
if missing_elements:
    raise ValueError(f"Error: 以下元素不存在于基准列表中: {missing_elements}")

for type_single in os.listdir(motif_dir):
    if type_single not in list(group['dir']):
        continue
    kr = os.path.join(motif_dir, type_single, "knownResults.txt")
    if not os.path.exists(kr):
        print(f'knownResults.txt not in {type_single}')
        continue
    df_temp = pd.read_table(kr)
    df_temp['Ratio'] = df_temp["% of Target Sequences with Motif"].map(
        lambda d: float(d[:-1])
    ) / df_temp["% of Background Sequences with Motif"].map(
        lambda d: float(d[:-1])
    )
    df_temp = df_temp[["Motif Name", "P-value", "Ratio"]]
    df_temp["Log P-value"] = df_temp["P-value"].map(lambda d: np.log10(d))
    for _, se in df_temp.iterrows():
        log_p_value_now = float(se["Log P-value"])
        r = float(se["Ratio"])
        m_n = se["Motif Name"].split("/", 1)[0]
        n_t = m_n.split("(", 1)
        t_t = name_handle(n_t[1]) if len(n_t) == 2 else ""
        s_t = type_single.rsplit("_", 1)
        data.append({
            #"sample": s_t[0].split('_')[0],
            #"cluster": s_t[0].split('_')[1],
            "dir": type_single,
            "motif": n_t[0],
            "type": t_t if len(t_t) > 0 else "Unknown",
            "Log P-value": log_p_value_now,
            "motif_name": m_n,
            "-log10(p-value)": -log_p_value_min if log_p_value_now < log_p_value_min else -log_p_value_now,
            "Ratio": r
        })
    del df_temp

if len(data) == 0:
    print("Error: No motif result found")
    sys.exit(1)

if not os.path.exists(save_dir):
    os.makedirs(save_dir)
df_data = pd.DataFrame(data)
df_data.sort_values(["dir", "type", "motif"], inplace=True)
df_data.to_csv("{}/homer_Allmotifs_result.txt".format(save_dir), sep="\t", index=False)

df_data = pd.merge(df_data, group, on='dir')
df_data['sample_cluster'] = df_data['sample']+'_'+df_data['cluster']

# 获得显著的motif，去掉不感兴趣的motif
y_need = set(df_data[df_data["Log P-value"] < log_p_value]["motif"].unique())
if choose_motifs != None:
    with open(choose_motifs, "r") as file:
    # 使用 readlines() 逐行读取并去除换行符
        chooseMotifs = [line.strip() for line in file.readlines()]
    y_need = set(chooseMotifs) & y_need
df_data = df_data[df_data["motif"].isin(y_need)]
df_data = df_data[df_data["Ratio"] != np.inf]

df_data_pivot = df_data.pivot_table(index="motif_name", columns="sample_cluster", values="-log10(p-value)")
df_data_pivot.fillna(0, inplace=True)

if len(df_data) == 0:
    print("motif no data, by filter")

# 为0的，后面画图用
df_mask = df_data_pivot.copy()
for i in df_mask.columns:
    df_mask[i] = df_mask[i].map(lambda d: d == 0)

group_y = dict()
group_x = dict()
for x, df_temp in df_data.groupby("sample"):
    group_x[x] = list(df_temp["sample_cluster"].unique())
    group_x[x].sort()
    del df_temp
for y, df_temp in df_data.groupby("type"):
    group_y[y] = list(df_temp["motif_name"].unique())
    group_y[y].sort()
    del df_temp
group_x_keys = list(group_x.keys())
group_x_keys.sort(reverse=True)
group_x_keys = tuple(group_x_keys)
group_y_keys = list(group_y.keys())
group_y_keys.sort()
group_y_keys = tuple(group_y_keys)

fig_size = (
    len(df_data['dir'].unique())*1+1,
    len(df_data['motif_name'].unique())*1+1
)
cmap = mpl.colors.LinearSegmentedColormap.from_list("www", ("white","#F2FAB0","#CE0013","#5E191A"))
norm = mpl.colors.Normalize(
    vmin=0,
    vmax=abs(log_p_value_min)
)

fig = plt.figure(figsize=fig_size)
grid = plt.GridSpec(
    len(df_data['motif_name'].unique())+3,
    len(df_data['sample_cluster'].unique()),
    hspace=0.2,
    wspace=0.2,
    figure=fig,
)

ax_list = list()
grid_min_y = 0
for y in group_y_keys:
    grid_min_x = 0
    grid_min_y_temp = grid_min_y + len(group_y[y])
    ax_list_temp = list()
    for i, x in enumerate(group_x_keys):
        grid_min_x_temp = grid_min_x + len(group_x[x])
        grid_temp = grid[grid_min_y:grid_min_y_temp, grid_min_x:grid_min_x_temp]
        if i == 0:
            if len(ax_list) == 0:
                ax_t = fig.add_subplot(
                    grid_temp,
                    xticklabels=[],
                    yticklabels=[],
                )
            else:
                ax_t = fig.add_subplot(
                    grid_temp,
                    xticklabels=[],
                    yticklabels=[],
                )
        else:
            if len(ax_list) == 0:
                ax_t = fig.add_subplot(
                    grid_temp,
                    # sharex=ax_list_temp[0],
                    xticklabels=[],
                    yticklabels=[],
                )
            else:
                ax_t = fig.add_subplot(
                    grid_temp,
                    # sharex=ax_list_temp[0],
                    xticklabels=[],
                    yticklabels=[],
                )
        ax_list_temp.append(ax_t)
        grid_min_x = grid_min_x_temp
    grid_min_y = grid_min_y_temp
    ax_list.append(ax_list_temp)

ax_list = list()
grid_min_y = 0
for y in group_y_keys:
    grid_min_x = 0
    grid_min_y_temp = grid_min_y + len(group_y[y])
    ax_list_temp = list()
    for i, x in enumerate(group_x_keys):
        grid_min_x_temp = grid_min_x + len(group_x[x])
        grid_temp = grid[grid_min_y:grid_min_y_temp, grid_min_x:grid_min_x_temp]
        if i == 0:
            if len(ax_list) == 0:
                ax_t = fig.add_subplot(
                    grid_temp,
                    xticklabels=[],
                    yticklabels=[],
                )
            else:
                ax_t = fig.add_subplot(
                    grid_temp,
                    xticklabels=[],
                    yticklabels=[],
                )
        else:
            if len(ax_list) == 0:
                ax_t = fig.add_subplot(
                    grid_temp,
                    # sharex=ax_list_temp[0],
                    xticklabels=[],
                    yticklabels=[],
                )
            else:
                ax_t = fig.add_subplot(
                    grid_temp,
                    # sharex=ax_list_temp[0],
                    xticklabels=[],
                    yticklabels=[],
                )
        ax_list_temp.append(ax_t)
        grid_min_x = grid_min_x_temp
    grid_min_y = grid_min_y_temp
    ax_list.append(ax_list_temp)

# g = max(int(grid_min_x * 0.2), 1)
# color_ax = fig.add_subplot(grid[grid_min_y+2:, g:-g], xticklabels=[], yticklabels=[])
if len(df_data['sample_cluster'].unique()) <= 2:
    g = 1  # 固定 g 的值
    color_ax = fig.add_subplot(grid[grid_min_y + 2:, :], xticklabels=[], yticklabels=[])
else:
    g = max(int(grid_min_x * 0.2), 1)
    color_ax = fig.add_subplot(grid[grid_min_y + 2:, g:-g], xticklabels=[], yticklabels=[])

mask_point = [[],[]]
log_p_value_point_temp = -log_p_value_point
for _, se in df_data.iterrows():
    # print(se["Log P-value"], se["Ratio"])
    if se["-log10(p-value)"] < log_p_value_point_temp:
        continue
    if se["Ratio"] < ratio_min:
        continue
    mask_point[0].append(se["sample_cluster"])
    mask_point[1].append(se["motif_name"])

def get_lable(ss, s, tick_index, label_index, s_n = 1):
    global name_handle
    tick_labels = list()
    label_name = set()
    for i in ss:
        if s in i:
            t = i.rsplit(s,s_n)
            tick_labels.append(t[tick_index])
            label_name.add(name_handle(t[label_index]))
        else:
            tick_labels.append(i)
    return tick_labels, label_name.pop() if len(label_name) == 1 else ""


def ax_single(y, x, y_max, x_max):
    global df_data_pivot, df_mask, cmap, norm, ax_list, group_x, group_y, group_x_keys, group_y_keys, mask_point, name_handle
    main_ax = ax_list[y][x]
    x_l = group_x_keys[x]
    y_l = group_y_keys[y]

    df_data_pivot_now = df_data_pivot.loc[group_y[y_l], group_x[x_l]]
    sns.heatmap(
        df_data_pivot_now,
        mask=df_mask.loc[list(df_data_pivot_now.index), list(df_data_pivot_now.columns)].to_numpy(),
        norm=norm,
        cmap=cmap,
        ax=main_ax,
        cbar=None,
        zorder=100,
        linecolor='#F0F0F0',
        linewidth=0.01,
        linestyle="--"
    )
    if (y_max - y) == 1:
        xtick_labels, xlabel_name = get_lable(df_data_pivot_now.columns, "_", 1, 0)
        main_ax.set_xticklabels(xtick_labels, rotation=0)
    else:
        xlabel_name = ""
        main_ax.xaxis.set_ticks([])
        main_ax.xaxis.set_major_locator(mpt.NullLocator())

    if x == 0:
        ytick_labels, ylabel_name = get_lable(df_data_pivot_now.index, "(", 0, 1)
        if x_max == 1:
            ylabel_name = name_handle(ylabel_name).replace(",", "\n")
            main_ax.yaxis.set_label_position("right")
            main_ax.set_ylabel(ylabel_name, rotation=0, va="center", ha="left")
        else:
            ylabel_name = ""
            main_ax.set_ylabel(ylabel_name, rotation=0)
        main_ax.set_yticklabels(ytick_labels, rotation=0)
    elif (x_max - x) == 1:
        ytick_labels, ylabel_name = get_lable(df_data_pivot_now.index, "(", 0, 1)
        ylabel_name = name_handle(ylabel_name).replace(",", "\n")
        main_ax.yaxis.set_ticks([])
        main_ax.yaxis.set_label_position("right")
        main_ax.set_ylabel(ylabel_name, rotation=0, va="center", ha="left")
    else:
        main_ax.yaxis.set_ticks([])
        ylabel_name = ""
        main_ax.set_ylabel(ylabel_name, rotation=0)

    # main_ax.grid(True, ls="--", lw=0.2, color='#F0F0F0', zorder=0)
    # main_ax.grid(True, ls="--", lw=0.5, color='black', zorder=0)

    x_half = 0
    x_s_start = df_data_pivot_now.columns[0].rsplit("_",1)[1][:2]
    for i in range(1, len(df_data_pivot_now.columns)):
        if df_data_pivot_now.columns[i].rsplit("_",1)[1][:2] == x_s_start:
            continue
        x_half = i
        break
    if (x_half > 0) and (x_half < (len(df_data_pivot_now.columns) - 1)):
        main_ax.plot(
            [x_half, x_half],
            [0, len(df_data_pivot_now.index)],
            ls="--",
            lw=0.8,
            c="white",
            zorder=110
        )

    y_dict = dict(map(lambda d: (d[1], d[0]), enumerate(df_data_pivot_now.index)))
    x_dict = dict(map(lambda d: (d[1], d[0]), enumerate(df_data_pivot_now.columns)))

    xxx = list()
    yyy = list()
    for i in range(len(mask_point[0])):
        if (mask_point[0][i] not in x_dict) or (mask_point[1][i] not in y_dict):
            continue
        xxx.append(x_dict[mask_point[0][i]]+0.5)
        yyy.append(y_dict[mask_point[1][i]]+0.5)
    main_ax.scatter(
        xxx,
        yyy,
        marker="*",
        c="black",
        zorder=110
    )

    main_ax.set_xlabel(xlabel_name)
    main_ax.yaxis.set_tick_params(left=False)
    main_ax.xaxis.set_tick_params(bottom=False)

    main_ax.spines[:].set_visible(True)
    main_ax.spines[:].set_zorder(110)

###
for y, x in product(range(len(group_y_keys)), range(len(group_x_keys))):
    ax_single(y, x, len(group_y_keys), len(group_x_keys))

def ax_single_corlor():
    global cmap, norm, log_p_value_point, log_p_value_min, color_ax

    plt.colorbar(
        mpm.ScalarMappable(
            norm=norm,
            cmap=cmap
        ),
        cax=color_ax,
        orientation='horizontal',
        label="-log10(p-value)",
    )
    color_ax.set_xticks(
        [0, -log_p_value_min],
        ["0", str(-log_p_value_min)]
    )

ax_single_corlor()

save_dir_img = os.path.join(save_dir,"img")
if not os.path.exists(save_dir_img):
    os.makedirs(save_dir_img)
plt.savefig(os.path.join(
    save_dir_img,
    "{}.pdf".format(
        os.path.basename(motif_dir)
    )
), dpi=500, bbox_inches='tight')
plt.close()

# return
save_dir_motif = os.path.join(save_dir, "motif")
if not os.path.exists(save_dir_motif):
    os.makedirs(save_dir_motif)
for y_i, y in enumerate(group_y_keys):
    ax_list = [[] for _ in range(y_i+1)]
    fig_size = (
        len(df_data['dir'].unique())*1+1,
        len(group_y[y])*1+1
    )
    fig = plt.figure(figsize=fig_size)
    grid = plt.GridSpec(
        len(group_y[y])+4,
        len(df_data['dir'].unique()),
        hspace=0.2,
        wspace=0.2,
        figure=fig,
    )
    grid_min_x = 0
    for i, x in enumerate(group_x_keys):
        grid_min_x_temp = grid_min_x + len(group_x[x])
        grid_temp = grid[:-4, grid_min_x:grid_min_x_temp]
        if i == 0:
            ax_t = fig.add_subplot(
                grid_temp,
                xticklabels=[],
                yticklabels=[],
            )
        else:
            ax_t = fig.add_subplot(
                grid_temp,
                sharex=ax_list_temp[0],
                yticklabels=[],
            )
        grid_min_x = grid_min_x_temp
        ax_list[-1].append(ax_t)
    if len(df_data['sample_cluster'].unique()) <= 2:
        g = 1  # 固定 g 的值
        color_ax = fig.add_subplot(grid[-1:, :], xticklabels=[], yticklabels=[])
    else:
        g = max(int(grid_min_x * 0.2), 1)
        color_ax = fig.add_subplot(grid[-1:, g:-g], xticklabels=[], yticklabels=[])
    #color_ax = fig.add_subplot(grid[-2:, g:-g], xticklabels=[], yticklabels=[])
    for x in range(len(group_x_keys)):
        ax_single(y_i, x, y_i+1, len(group_x_keys))
    ax_single_corlor()
    # plt.show()

    plt.savefig(os.path.join(
        save_dir_motif,
        "{}_{}.pdf".format(
            os.path.basename(motif_dir),
            y
        )
    ), dpi=500, bbox_inches='tight')
    plt.close()

print('运行完成，可以查看结果了！')
