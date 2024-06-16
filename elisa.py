import os
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import xlrd
import csv
from itertools import compress
import shutil


# Functions


def read_plate(sheet):
    """
        Given a spreadsheet of a plate, read it and extract the data

        :param sheet: The workspace file path

        :return: patient_id: Patient ID
        :return: sample_id: Sample ID
        :return: cytokine: Analyzed cytokine name
        :return: plate_values: dictionary to hold target labels and the respective ODs
    """
    # Get patient ID and cytokine tested
    patient_id = int(sheet.cell_value(0, 1))
    sample_id = int(sheet.cell_value(1, 1))
    cytokine = sheet.cell_value(2, 1)

    # Plate layout dictionary
    plate_layout = {}

    # How many rows the plate layout table is ahead of the start
    layout_row_offset = 13

    # How many columns the plate layout table (and the other tables) is ahead of the start
    col_offset = 2

    # Fill the dictionary with the row and column numbers of cells containing each label
    for row in range(0, 8):
        for col in range(0, 12):
            label = sheet.cell_value(row + layout_row_offset, col + col_offset)
            if label not in plate_layout:
                plate_layout[label] = [[row, col]]
            else:
                plate_layout[label].append([row, col])

    # How many rows the values table is ahead of the start
    values_row_offset = 4

    # How many rows the dilution factor table is ahead of the start
    dilution_row_offset = 22

    # Dictionary to hold plate labels and respective values
    plate_values = {}

    # Fill values dictionary with label as key and respective average value times the respective dilution factor
    for label, indices in plate_layout.items():
        sum = 0
        for cell in indices:
            sum += sheet.cell_value(cell[0] + values_row_offset, cell[1] + col_offset)
        value = sum / len(indices)
        dilution_factor = sheet.cell_value(cell[0] + dilution_row_offset, cell[1] + col_offset)
        plate_values[label] = value * dilution_factor

    # If the medium exists as control in the plate, subtract it from all other mappings (except Standards)
    if 'medium' in plate_layout:
        for k, v in plate_values.items():
            if not k.startswith('STD') and k != 'medium':
                plate_values[k] = v - plate_values['medium']

    return patient_id, sample_id, cytokine, plate_values


# Curve fit function
def log4pl(x, A, B, C, D):
    """
        4 parameter logistic regression function

        :param x: value
        :param A: minimum value of OD
        :param B: slope of the curve at point C
        :param C: point of inflection
        :param D: maximum value of OD
        :return: concentration value
    """
    return ((A - D) / (1.0 + ((x / C) ** B))) + D


def to_pct(value, max):
    """
        Convert the given value to percentage, where max is the reference for 100%

        :param value: value
        :param max: maximum value
        :return: value in percentage
    """
    return value / max * 100


def calculate_change(medium, peptide):
    return round(float(peptide * 100 / medium), 3)


# ---------------------------------------------------------------------------------
# Execution starts here

if os.path.exists(f'{os.getcwd()}/Patients'):
    shutil.rmtree(f'{os.getcwd()}/Patients', ignore_errors=True)

controls = ['PHA', 'OKT3', 'medium']
antigens_controls = ['CMV', 'EBNA', 'M1 (mix)', 'ESAT6', 'Haemagluttinin']

base_dir = "Data/DATA_Raw_files/ELISA"
analysis_type = ''

if not os.path.exists('Patients'):
    os.mkdir('Patients')

for dir in os.listdir(base_dir):
    if str(dir) == 'WBA':
        analysis_type = 'WBA'
    else:
        analysis_type = 'TIL'

    for file in os.listdir(f'{base_dir}/{dir}'):
        if file.endswith('xls'):
            wb = xlrd.open_workbook(f'{base_dir}/{dir}/{file}')
            for i in range(len(wb.sheet_names())):
                sheet = wb.sheet_by_index(i)
                # Get the plate values
                patient_id, sample_id, cytokine, plate_values = read_plate(sheet)

                print(f'Analysis Type: {analysis_type} | Patient {patient_id} | Sample {sample_id}')

                # Get standards
                plate_values_std = {k: v for k, v in plate_values.items() if k.startswith('STD')}

                # Filter out targets that did not react
                targets = {k: v for k, v in plate_values.items() if v > 0 and not k.startswith('STD') and k != 'medium'}

                print(f'Targets Used: {list(plate_values.keys())}')
                print(f'Targets That Reacted: {list(targets.keys())}\n')

                filepath = f'Patients/Patient {patient_id}/{analysis_type}/Sample {sample_id}/Cytokine {cytokine}'

                if not os.path.exists(filepath):
                    os.makedirs(filepath)

                # Build standard curve
                x = [1000, 500, 250, 125, 62.5, 31.25, 15.625, 7.8125]  # Concentration
                y = list(plate_values_std.values())  # Optical Density

                # X values to draw the curve
                x_min = min(x)
                x_max = max(x)
                step = 0.01
                x_new = np.arange(x_min, x_max, step)

                params, _ = curve_fit(log4pl, x, y)
                A, B, C, D = params[0], params[1], params[2], params[3]
                print("4PL parameters: A = " + str(round(A, 4)) + ", B = " + str(round(B, 4)) +
                      ", C = " + str(round(C, 4)) + ", D = " + str(round(D, 4)))

                # Y values to draw the curve
                yfit1_new = ((A - D) / (1.0 + ((x_new / C) ** B))) + D

                # Plot the curve
                plt.plot(x, y, 'r+', label="y-original")
                plt.plot(x_new, yfit1_new, label="y=((A-D)/(1.0+((x/C)^B))) + D")
                plt.xlabel('OD')
                plt.ylabel('Standards (log)')
                plt.xscale("log")
                plt.legend(loc='best', fancybox=True, shadow=True)
                plt.grid(True)
                plt.savefig(f'{filepath}/{analysis_type}_standard_curve.png')
                plt.clf()

                absError = log4pl(x, A, B, C, D) - y
                SE = np.square(absError)  # squared errors
                MSE = np.mean(SE)  # mean squared errors
                RMSE = np.sqrt(MSE)  # Root Mean Squared Error, RMSE
                Rsquared = 1.0 - (np.var(absError) / np.var(y))
                print('RMSE:', RMSE)
                print('R-squared:', Rsquared)
                print()

                # Plot bar chart
                concentrations = {k: round(log4pl(v, A, B, C, D), 4) for k, v in targets.items()}
                data_min = min(concentrations.values())
                data_max = max(concentrations.values())

                fig, ax = plt.subplots(figsize=(16, 9))

                labels = list(concentrations.keys())
                values = list(concentrations.values())

                plt.xlabel(f'{cytokine} pg/mL', fontsize=15)

                controls_mask = [bool(item in controls) for item in labels]
                antigens_mask = [bool(item in antigens_controls) for item in labels]
                none_mask = [not (a or b) for a, b in zip(controls_mask, antigens_mask)]

                ax.barh(list(compress(labels, controls_mask)), list(compress(values, controls_mask)), color='blue')
                ax.barh(list(compress(labels, antigens_mask)), list(compress(values, antigens_mask)), color='red')
                ax.barh(list(compress(labels, none_mask)), list(compress(values, none_mask)), color='green')

                plt.yticks(fontsize=12)

                # Add padding between axes and labels
                ax.xaxis.set_tick_params(pad=5)
                ax.yaxis.set_tick_params(pad=10)

                # Add x, y gridlines
                ax.grid(color='grey',
                        linestyle='-.', linewidth=0.5,
                        alpha=0.2)

                # Add annotation to bars
                for i in ax.patches:
                    plt.text(i.get_width(), i.get_y() + 0.35,
                             f' {str(round((i.get_width()), 4))}',
                             fontsize=10, fontweight='bold',
                             color='grey')

                plt.title(f'Patient {patient_id} - Sample {sample_id}', fontsize=25)

                red_patch = mpatches.Patch(color='red', label='Viral Antigens')
                blue_patch = mpatches.Patch(color='blue', label='Controls')
                green_patch = mpatches.Patch(color='green', label='Peptides')

                plt.legend(handles=[red_patch, blue_patch, green_patch])

                # Show Plot
                plt.subplots_adjust(left=0.2)
                plt.xlim(data_min - data_min/100, data_max + data_max/50)

                plt.savefig(f'{filepath}/{analysis_type}_concentrations.png')
                plt.clf()
                plt.subplots_adjust(left=0.1)

                # Get ODs in order to calculate relative percentages
                sorted_ods = sorted(targets.items(), key=lambda x: x[1])
                sorted_ods = dict(sorted_ods)

                od_values = list(sorted_ods.values())
                max_value = od_values[len(od_values) - 1]

                sorted_ods = {k: round(to_pct(v, max_value), 2) for k, v in sorted_ods.items()}

                peptide_reactions = {k: calculate_change(plate_values['medium'], v) for k, v in targets.items()}

                print(peptide_reactions)

                # Write Peptide reactions data to sample file
                with open(os.path.join(filepath, f'peptide_reactions.csv'), 'w', newline='\n') as f:
                    w = csv.writer(f)
                    w.writerow(peptide_reactions.keys())
                    w.writerow(peptide_reactions.values())

                # Write Peptide reactions data to global file
                with open(f'{os.getcwd()}/Patients/{analysis_type}_global_peptide_reactions.csv', 'a', newline='') as f:
                    w = csv.writer(f)
                    w.writerow(peptide_reactions.keys())
                    w.writerow(peptide_reactions.values())

                # Remove controls
                sorted_ods = {k: v for k, v in sorted_ods.items() if k not in controls and k not in antigens_controls}
                concentrations = {k: v for k, v in concentrations.items() if k not in controls and k not in antigens_controls}

                print(f'Sorted ODs: {sorted_ods}\n')
                print(f'Concentrations: {concentrations}')

                # Write OD data to sample file
                with open(os.path.join(filepath, f'{analysis_type}_ODs.csv'), 'w', newline='\n') as f:
                    w = csv.writer(f)
                    w.writerow(sorted_ods.keys())
                    w.writerow(sorted_ods.values())

                # Write OD data to global file
                with open(f'{os.getcwd()}/Patients/{analysis_type}_global_ODs.csv', 'a', newline='') as f:
                    w = csv.writer(f)
                    w.writerow(sorted_ods.keys())
                    w.writerow(sorted_ods.values())

                # Write concentrations data to sample file
                with open(os.path.join(filepath, f'{analysis_type}_Concentrations.csv'), 'w', newline='\n') as f:
                    w = csv.writer(f)
                    w.writerow(labels)
                    w.writerow(values)

                # Write concentrations data to global file
                with open(f'{os.getcwd()}/Patients/{analysis_type}_global_concentrations.csv', 'a', newline='') as f:
                    w = csv.writer(f)
                    w.writerow(labels)
                    w.writerow(values)

