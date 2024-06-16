from scipy.stats import linregress
import numpy as np
from global_func import get_data, get_gate_data
import matplotlib.pyplot as plt
from scipy.stats import t
import os
import csv


# Functions

def get_avg(peptide, data):
    """
        Get the average mfi for a peptide in the data provided
    :param peptide: the peptide string
    :param data: the data to be checked
    :return: the average mfi of the peptide in the data list
    """
    res = []
    for dict_item in data:
        if peptide in dict_item:
            res.append(float(dict_item[peptide]))
    return round(sum(res) / len(res), 3)


def perform_linear_regression(data_min, data_max, elisa_data, analysis, gate, elisa_type):
    """
        Performs linear regression with the given data

    :param data_min: the minimum value of the gate data
    :param data_max: the maximum value of the gate data
    :param elisa_data: the elisa data
    :param analysis: analysis string (PBMCs or TILs)
    :param gate: the gate (CD4, CD8, GD)
    :param elisa_type: the ELISA analysis type (TIL or WBA)
    :return: None
    """
    reg_array = np.random.uniform(low=float(data_min), high=float(data_max), size=(len(elisa_data),))
    res = linregress(reg_array, elisa_data)

    print(f'R-squared: {res.rvalue**2:.6f}')

    tinv = lambda p, df: abs(t.ppf(p / 2, df))
    ts = tinv(0.05, len(reg_array) - 2)

    print(f"slope (95%): {res.slope:.6f} +/- {ts * res.stderr:.6f}")
    print(f"intercept (95%): {res.intercept:.6f}"
          f" +/- {ts * res.intercept_stderr:.6f}")

    extracted_data = {'R-Squared': round(res.rvalue**2, 6),
                      'Slope': f'{round(res.slope, 6)} +/- {round(ts + res.stderr, 6)}',
                      'Intercept': f'{round(res.intercept, 6)} +/- {round(ts * res.intercept_stderr, 6)}'}

    filepath = f'LinearRegression/{gate}_IFN-g_{analysis}_vs_ELISA_{elisa_type}'
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    # Write data to sample file
    with open(f'{filepath}/data.csv', 'w', newline='\n') as f:
        w = csv.writer(f)
        w.writerow(extracted_data.keys())
        w.writerow(extracted_data.values())

    plt.plot(reg_array, elisa_data, 'o', label='original data')
    plt.plot(reg_array, res.intercept + res.slope * reg_array, 'r', label='fitted line')

    plt.title(f'{gate} IFN-g {analysis} VS ELISA {elisa_type}', fontsize=15)

    plt.xticks([])
    plt.yticks([])

    plt.legend()

    plt.savefig(f'{filepath}/graph.png')
    plt.clf()

    print()


# ---------------------------------------------------------------------------------
# Execution starts here

# If results directorydoes not exist, create it
if not os.path.exists('LinearRegression'):
    os.makedirs('LinearRegression')

# get elisa data
elisa_wba = get_data('ELISA', 'Reactions', 'WBA', '')
elisa_til = get_data('ELISA', 'Reactions', 'TIL', '')

print(f'Elisa WBA: {elisa_wba}')
print(f'Elisa TIL: {elisa_til}')

# Rebuild elisa data into a dictionary to take only the average of the peptides it holds
final_elisa_wba = {}
for dict_item in elisa_wba:
    for k in dict_item.keys():
        if k not in final_elisa_wba:
            final_elisa_wba[k] = get_avg(k, elisa_wba)

final_elisa_til = {}
for dict_item in elisa_til:
    for k in dict_item.keys():
        if k not in final_elisa_til:
            final_elisa_til[k] = get_avg(k, elisa_til)

print(f'Elisa WBA Values: {list(final_elisa_wba.values())}')
print(f'Elisa TIL Values: {list(final_elisa_til.values())}')
print()

# PBMCs gate percentages vs ELISA WBA
# CD8
cd8_pbmcs_data = [float(list(entry.values())[0]) for entry in get_gate_data('cd8', 'FLOW', 'ICS', 'PBMCs', 'gate_pct')]
perform_linear_regression(min(cd8_pbmcs_data), max(cd8_pbmcs_data), list(final_elisa_wba.values()), 'PBMCs', 'CD8', 'WBA')

# CD4
cd4_pbmcs_data = [list(entry.values())[0] for entry in get_gate_data('cd4', 'FLOW', 'ICS', 'PBMCs', 'gate_pct')]
perform_linear_regression(min(cd4_pbmcs_data), max(cd4_pbmcs_data), list(final_elisa_wba.values()), 'PBMCs', 'CD4', 'WBA')
# GD
gd_pbmcs_data = [list(entry.values())[0] for entry in get_gate_data('gd', 'FLOW', 'ICS', 'PBMCs', 'gate_pct')]
perform_linear_regression(min(gd_pbmcs_data), max(gd_pbmcs_data), list(final_elisa_wba.values()), 'PBMCs', 'GD', 'WBA')

# TILs gate percentages vs ELISA TIL
# CD8
cd8_tils_data = [list(entry.values())[0] for entry in get_gate_data('cd8', 'FLOW', 'ICS', 'TILs', 'gate_pct')]
perform_linear_regression(min(cd8_tils_data), max(cd8_tils_data), list(final_elisa_til.values()), 'TILs', 'CD8', 'TIL')

# CD4
cd4_tils_data = [list(entry.values())[0] for entry in get_gate_data('cd4', 'FLOW', 'ICS', 'TILs', 'gate_pct')]
perform_linear_regression(min(cd4_tils_data), max(cd4_tils_data), list(final_elisa_til.values()), 'TILs', 'CD4', 'TIL')

# GD
gd_tils_data = [list(entry.values())[0] for entry in get_gate_data('gd', 'FLOW', 'ICS', 'TILs', 'gate_pct')]
perform_linear_regression(min(gd_tils_data), max(gd_tils_data), list(final_elisa_til.values()), 'TILs', 'GD', 'TIL')

