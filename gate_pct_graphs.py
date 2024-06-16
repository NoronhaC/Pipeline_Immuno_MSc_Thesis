import math

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

import global_func
import statistics


def get_comp_data(gate_list, group):
    avg = {}
    print(gate_list)
    print(group)
    get_group_avg(avg, gate_list, group)
    print(avg)
    return avg


def get_group_avg(avg, gate_list, group):
    for gate in gate_list:
        data_til = global_func.get_gate_data(gate, 'FLOW', group, 'TILs', 'gate_pct')
        data_pbmc = global_func.get_gate_data(gate, 'FLOW', group, 'PBMCs', 'gate_pct')

        avg[gate] = {'TILs': {'mean': 0, 'std': 0}, 'PBMCs': {'mean': 0, 'std': 0}}

        if data_til:
            data_values = [float(value) for dict_item in data_til for value in list(dict_item.values())]
            avg[gate]['TILs']['mean'] = round(statistics.mean(data_values), 3)
            avg[gate]['TILs']['std'] = round(statistics.stdev(data_values), 3)
        if data_pbmc:
            data_values = [float(value) for dict_item in data_pbmc for value in list(dict_item.values())]
            avg[gate]['PBMCs']['mean'] = round(statistics.mean(data_values), 3)
            avg[gate]['PBMCs']['std'] = round(statistics.stdev(data_values), 3)


ics_gate_list = ['CD4+', 'CD8+']
tcell_gate_list = ['precursor', 'TCM', 'TEM', 'TEFF', 'CD103+', 'CD95+', 'CD39+ CD69+']
thoming_gate_list = ['Th1', 'Th2', 'Th17', 'Tscm']

ics_data = get_comp_data(ics_gate_list, 'ICS')
tcell_data = get_comp_data(tcell_gate_list, 'Tcell')
thoming_data = get_comp_data(thoming_gate_list, 'Thoming')

final_data = {'ICS': ics_data, 'Tcell': tcell_data, 'Thoming': thoming_data}

for key, dataset in final_data.items():
    print(dataset)

    _, ax = plt.subplots(figsize=(16, 9))

    if key == 'ICS':
        plt.title(f'{key} IFN-g', fontsize=30)
    else:
        plt.title(f'{key}', fontsize=30)

    legend = ['PBMCs', 'TILs']
    custom_lines = [Line2D([0], [0], color="cyan", lw=6),
                    Line2D([0], [0], color="purple", lw=6)]

    ax.legend(custom_lines, legend, loc='best', fontsize=15)

    gates = list(dataset.keys())
    if key == 'Tcell':
        gates = ['CD8+ | ' + gate for gate in gates]
    elif key == 'Thoming':
        gates = ['CD4+ | ' + gate if gate.lower() != 'tscm' else 'CD8+ | ' + gate for gate in gates]

    x = np.arange(len(dataset))
    width = 0.4

    pbmcs_avg = []
    pbmcs_std = []
    tils_avg = []
    tils_std = []
    for values in list(dataset.values()):
        pbmcs_avg.append(values['PBMCs']['mean'])
        pbmcs_std.append(values['PBMCs']['std'])
        tils_avg.append(values['TILs']['mean'])
        tils_std.append(values['TILs']['std'])

    ax.bar(x - 0.2, pbmcs_avg, width, color='cyan')
    ax.bar(x + 0.2, tils_avg, width, color='purple')

    all_std = pbmcs_std + tils_std
    iterator = iter(all_std)

    for p in ax.patches:
        if p.get_height() > 0:
            rect_x = p.get_x()  # get the bottom left x corner of the bar
            rect_y = p.get_y()
            w = p.get_width()  # get width of bar
            h = p.get_height()  # get height of bar
            rect_std = next(iterator)
            min_y = h - rect_std  # use h to get min from dict z
            max_y = h + rect_std  # use h to get max from dict z
            plt.vlines(rect_x + w / 2, min_y, max_y, color='k')  # draw a vertical line
            plt.hlines(max_y, (rect_x + w / 2) - w / 4, (rect_x + w / 2) + w / 4, color='k')
            plt.hlines(min_y, (rect_x + w / 2) - w / 4, (rect_x + w / 2) + w / 4, color='k')

    plt.xlabel("Gates", fontsize=20)
    plt.ylabel("Gate Percentages", fontsize=20)

    plt.xticks(x, gates, rotation=30, ha="center", fontsize=20)
    plt.yticks(fontsize=20)
    plt.subplots_adjust(bottom=0.25)

    plt.savefig(f'Samples/{key}_gate_percentages.png')
    plt.clf()
