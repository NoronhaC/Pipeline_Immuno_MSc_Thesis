import scipy.stats as stats
from global_func import get_data
from tkinter import simpledialog
import os


def get_average(key, th_list):
    """
        Return the average of the gate mfi for the given Thoming list

        :param key: the gate to be checked
        :param th_list: th1 or th2 list
        :return: the average of the gate specified in the th list specified
    """
    res = []
    for dict_item in th_list:
        dict_item = eval(dict_item)
        if key in dict_item:
            res.append(dict_item[key])
    return round(sum(res) / len(res), 3)


def get_values_for_key(key):
    """
        filter the averages dictionary given the key (Th1 or Th2)
    :param key: Th1 or Th2
    :return: Dictionary with PBMCs and TILs for the given key
    """
    return {k: v for k, v in averages.items() if k.startswith(key)}

# ---------------------------------------------------------------------------------
# Execution starts here


# Use tkinter to get user input for significance level (alpha) for p-value comparison
alpha = simpledialog.askfloat('T-Test Dependent Tests',
                              'Input alpha value for p-value comparisons:')

if alpha is None or alpha < 0:
    print('ERROR: Alpha value is negative or not defined.')
    exit()

input('Press ENTER to continue...')

th1_gates = ['CXCR3', 'CCR6']
th2_gates = ['CCR4', 'CCR6']
all_gates = list(set(th1_gates + th2_gates))

th1_pbmcs = [dict_item['Th1'] for dict_item in get_data('FLOW', 'Thoming', 'PBMCs', 'gate_mfi') if 'Th1' in dict_item]
th2_pbmcs = [dict_item['Th2'] for dict_item in get_data('FLOW', 'Thoming', 'PBMCs', 'gate_mfi') if 'Th2' in dict_item]

th1_tils = [dict_item['Th1'] for dict_item in get_data('FLOW', 'Thoming', 'TILs', 'gate_mfi') if 'Th1' in dict_item]
th2_tils = [dict_item['Th2'] for dict_item in get_data('FLOW', 'Thoming', 'TILs', 'gate_mfi') if 'Th2' in dict_item]

# Get the averages for all gates composing each Thoming gate for the respective analysis type
averages = {'Th1_PBMCs': [get_average('CXCR3', th1_pbmcs), get_average('CCR6', th1_pbmcs)],
            'Th2_PBMCs': [get_average('CCR4', th2_pbmcs), get_average('CCR6', th2_pbmcs)],
            'Th1_TILs': [get_average('CXCR3', th1_tils), get_average('CCR6', th1_tils)],
            'Th2_TILs': [get_average('CCR4', th2_tils), get_average('CCR6', th2_tils)]
            }

print(averages)

th1_pbmcs_norm = stats.anderson(averages['Th1_PBMCs'], dist='norm').fit_result.success
th1_tils_norm = stats.anderson(averages['Th1_TILs'], dist='norm').fit_result.success

th2_pbmcs_norm = stats.anderson(averages['Th2_PBMCs'], dist='norm').fit_result.success
th2_tils_norm = stats.anderson(averages['Th2_TILs'], dist='norm').fit_result.success

if th1_pbmcs_norm and th1_tils_norm:
    # Perform related T-test for Th1
    th1 = get_values_for_key('Th1')
    th1_t_stat, th1_p_value = stats.ttest_rel(th1['Th1_PBMCs'], th1['Th1_TILs'])
    print(f'Th1:\nT-stat: {th1_t_stat}\t\tP-value:{th1_p_value}')

    if th1_p_value <= alpha:
        print('Reject the null hypothesis: the averages of CXCR3 and CCR6 for Th1 PBMCs and Th1 TILs are different')
    else:
        print('Accept the null hypothesis: the averages of CXCR3 and CCR6 for Th1 PBMCs and Th1 TILs are similar')

else:
    th1_pbmcs_wilcoxon = stats.wilcoxon(averages['Th1_PBMCs'])
    print(th1_pbmcs_wilcoxon)
    th1_tils_wilcoxon = stats.wilcoxon(averages['Th1_TILs'])
    print(th1_tils_wilcoxon)


if th2_pbmcs_norm and th2_tils_norm:
    # Perform related T-test for Th2
    th2 = get_values_for_key('Th2')
    th2_t_stat, th2_p_value = stats.ttest_rel(th2['Th2_PBMCs'], th2['Th2_TILs'])
    print(f'Th2:\nT-stat: {th2_t_stat}\t\tP-value:{th2_p_value}')

    if th2_p_value <= alpha:
        print('Reject the null hypothesis: the averages of CCR4 and CCR6 for Th2 PBMCs and Th1 TILs are different')
    else:
        print('Accept the null hypothesis: the averages of CCR4 and CCR6 for Th2 PBMCs and Th1 TILs are similar')

else:
    th2_pbmcs_wilcoxon = stats.wilcoxon(averages['Th2_PBMCs'])
    print(th2_pbmcs_wilcoxon)
    th2_tils_wilcoxon = stats.wilcoxon(averages['Th2_TILs'])
    print(th2_tils_wilcoxon)