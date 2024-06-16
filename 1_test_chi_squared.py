import scipy.stats as stats
import os
from global_func import get_data, get_gate_data
from tkinter import simpledialog


# Functions


def print_comp_data(comp_str, list1, list2):
    """
            Print both lists of data

            :param comp_str: The string that describes the two data lists
            :param list1: the first list
            :param list2: the second list

            :return: None
    """
    print(comp_str)
    print()
    print(list1)
    print()
    print(list2)
    print()


def perform_chi2_test(list1, list2, comp_str):
    """
            Perform the chi-squared contingency test with the
            values on the two lists

            :param list1: the first list
            :param list2: the second list
            :param comp_str: The string that describes the two data lists

            :return: None
    """
    # Create the directory structure
    if not os.path.exists('Stats_Tests/Chi2'):
        os.makedirs('Stats_Tests/Chi2')

    written_ind = False
    written_dep = False

    for dict_item in list1:
        for other_dict_item in list2:
            for key1, value1 in dict_item.items():
                for key2, value2 in other_dict_item.items():
                    value1 = float(value1)
                    value2 = float(value2)
                    print(f'Chi-Squared test for {key1}:{value1} and {key2}:{value2}')
                    if value1 != 100 or value2 != 100:
                        values_1 = [value1, float(100 - value1)]
                        values_2 = [value2, float(100 - value2)]
                        chi2, pvalue, dof, expected_freq = stats.chi2_contingency([values_1, values_2])
                        print(f'Chi-Squared: {chi2}\t\tP-Value: {pvalue}\t\tDegrees of Freedom: {dof}\t\t'
                              f'Expected Frequency: {expected_freq}')

                        if pvalue <= alpha:
                            print(f'dependent: the relation between {key1} and {key2} is significant')
                            file = open(f'Stats_Tests/Chi2/results_chi2_dependent_{alpha}.csv', 'a')
                            if not written_dep:
                                file.write(comp_str + '\n')
                                written_dep = True

                            file.write(f'{key1} versus {key2}')
                            file.write('\n')
                            file.write(f'Chi2: {chi2}, P-Value: {pvalue}, DOF: {dof}, Expected_Freq: {expected_freq}')
                            file.write('\n')
                        else:
                            print(f'independent: {key1} and {key2} do not have a significant relation')
                            file = open(f'Stats_Tests/Chi2/results_chi2_independent_{alpha}.csv', 'a')
                            if not written_ind:
                                file.write(comp_str + '\n')
                                written_ind = True

                            file.write(f'{key1} versus {key2}')
                            file.write('\n')
                            file.write(f'Chi2: {chi2}, P-Value: {pvalue}, DOF: {dof}, Expected_Freq: {expected_freq}')
                            file.write('\n')
            break
    file.write('\n')
    file.close()


# ---------------------------------------------------------------------------------
# Execution starts here

# Use tkinter to get user input for significance level (alpha) for p-value comparison
alpha = simpledialog.askfloat('Chi-Squared Contingency Tests',
                              'Input alpha value for p-value comparisons:')

if alpha is None or alpha < 0:
    print('ERROR: Alpha value is negative or not defined.')
    exit()

input('Press ENTER to continue...')

# If the files to be created already exist, delete them
if os.path.exists(f'Stats_Tests/Chi2/results_chi2_dependent_{alpha}.csv'):
    os.remove(f'Stats_Tests/Chi2/results_chi2_dependent_{alpha}.csv')
if os.path.exists(f'Stats_Tests/Chi2/results_chi2_independent_{alpha}.csv'):
    os.remove(f'Stats_Tests/Chi2/results_chi2_independent_{alpha}.csv')

# gate % flow_Thoming_PBMCs                 versus  gate % flow_Thoming_TILs
flow_thoming_pbmcs_gate_pct = get_data('FLOW', 'Thoming', 'PBMCs', 'gate_pct')
flow_thoming_tils_gate_pct = get_data('FLOW', 'Thoming', 'TILs', 'gate_pct')
print_comp_data('GATE PCT THOMING PBMCS VS GATE PCT THOMING TILS', flow_thoming_pbmcs_gate_pct,
                flow_thoming_tils_gate_pct)
perform_chi2_test(flow_thoming_pbmcs_gate_pct, flow_thoming_tils_gate_pct, 'Thoming PBMCs VS Thoming TILs')

# gate % flow_Tcell_PBMCs                 versus  %(each peptide) ELISA_WBA
flow_tcell_pbmcs_gate_pct = get_data('FLOW', 'Tcell', 'PBMCs', 'gate_pct')
elisa_wba = get_data('ELISA', 'OD', 'WBA', '')
print_comp_data('GATE PCT TCELL PBMCS VS ELISA', flow_tcell_pbmcs_gate_pct, elisa_wba)
perform_chi2_test(flow_tcell_pbmcs_gate_pct, elisa_wba, 'Tcell PBMCs VS ELISA WBA')

# gate % flow_ICS_PBMCs                   versus  %(each peptide) ELISA_WBA
flow_ics_pbmcs_gate_pct = get_data('FLOW', 'ICS', 'PBMCs', 'gate_pct')
print_comp_data('GATE PCT ICS PBMCS VS ELISA', flow_ics_pbmcs_gate_pct, elisa_wba)
perform_chi2_test(flow_ics_pbmcs_gate_pct, elisa_wba, 'ICS PBMCs VS ELISA WBA')

# 	ICS_PBMCs_CD3_IFNgamma   COMPARE with ELISA_WBA
print('MFI CD3 ICS PBMCS VS ELISA')
print('No data.')
print()

# 	ICS_PBMCs_CD8_IFNgamma   COMPARE with ELISA_WBA
cd8_data = get_gate_data('cd8', 'FLOW', 'ICS', 'PBMCs', 'gate_pct')
print_comp_data('GATE PCT CD8 ICS PBMCS VS ELISA', cd8_data, elisa_wba)
perform_chi2_test(cd8_data, elisa_wba, 'ICS PBMCS CD8 INF-G VS ELISA WBA')

# 	ICS_PBMCs_CD4_IFNgamma   COMPARE with ELISA_WBA
cd4_data = get_gate_data('cd4', 'FLOW', 'ICS', 'PBMCs', 'gate_pct')
print_comp_data('GATE PCT CD4 ICS PBMCS VS ELISA', cd4_data, elisa_wba)
perform_chi2_test(cd4_data, elisa_wba, 'ICS PBMCS CD4 INF-G VS ELISA WBA')

# 	ICS_PBMCs_gd_IFNgamma    COMPARE with ELISA_WBA
gd_data = get_gate_data('gd', 'FLOW', 'ICS', 'PBMCs', 'gate_pct')
print_comp_data('GATE PCT GD ICS PBMCS VS ELISA', gd_data, elisa_wba)
perform_chi2_test(gd_data, elisa_wba, 'ICS PBMCS GD INF-G VS ELISA WBA')
