# NOTE: This script requires graphviz to be installed in order to run: https://graphviz.org/download/

import os
import flowkit as fk
import csv
import graphviz
import shutil


# Functions

def get_fluoro_labels(sample, events):
    """
        Extract fluorescent event labels from a sample

        :param sample: The sample instance
        :param events: The gate events
        :return: an array containing the names of the columns that correspond to the fluorescent labels
    """
    fluoro_labels = []
    for chan_idx in sample.fluoro_indices:
        fluoro_labels.append(sample.pnn_labels[chan_idx])

    events_labels = []
    for label in fluoro_labels:
        for other_label in list(events.columns):
            if label in other_label:
                events_labels.append(other_label)
    return events_labels


def get_wsp_files(dir):
    """
        Recursive function to get all the workspace files from a directory

        :param dir: The base directory
        :return: an array with the paths to all workspace files
    """
    wsp_files = []
    for item in os.listdir(dir):
        item_path = os.path.join(dir, item)
        # if the item is a directory, get all workspace files from that directory
        if os.path.isdir(item_path):
            aux_list = get_wsp_files(item_path)
            for file in aux_list:
                wsp_files.append(file)
        # else, if it is a file, add it to the list if it is a workspace file
        elif os.path.isfile(item_path):
            if item.endswith('.wsp'):
                wsp_files.append(item_path)
    return wsp_files


def extract_data(gate, gate_path, gate_pct, mfi_comp, parent, sample, sample_id, sample_results, wsp, labels_dict):
    """
        Refactored piece of code that extracts the MFI for required gates and gate percentages
        and appends them to the respective lists

        :param gate: The gate
        :param gate_path: A tuple that represents the path to the gate
        :param gate_pct: The gate percentages list
        :param mfi_comp: The MFI comparison list
        :param parent: The parent gate
        :param sample: The sample instance
        :param sample_id: The ID of the sample instance
        :param sample_results: The sample results dataframe
        :param wsp: the workspace file
        :param labels_dict: labels dictionary

        :return: None
    """

    # Get the index for the correct report entry
    idx = -1
    for item in list(sample_results.report.loc[sample_results.report['gate_path'] == gate_path].gate_name):
        idx += 1
        if gate == item.strip():
            break

    # retrieve the relative percent from the report that matches the gate being handled
    relative_pct = float(sample_results.report.loc[sample_results.report['gate_path'] == gate_path].iloc[[idx]]
                         .relative_percent)

    # add the value to the gate percentages dictionary
    gate_name = ''
    if pre == 'Thoming':
        if th1_gate(gate=gate, gate_path=gate_path, parent=parent):
            gate_name = gate + " Th1" if "th1" not in gate.lower() else gate
        elif th2_gate(gate, gate_path, parent):
            gate_name = gate + " Th2" if "th2" not in gate.lower() else gate
        elif th17_gate(gate, gate_path, parent):
            gate_name = gate + " Th17" if "th17" not in gate.lower() else gate
        elif th1star_gate(gate, gate_path, parent):
            gate_name = gate + " Th1*" if "th1*" not in gate.lower() or "th1+" not in gate.lower() else gate
        elif tscm_gate(gate, gate_path, parent):
            gate_name = gate + " Tscm" if "tscm" not in gate.lower() else gate
        gate_pct[f'{parent} | {gate_name}'] = round(relative_pct, 3)

    gate_pct[f'{parent} | {gate}'] = round(relative_pct, 3)

    # Now to calculate the MFI
    gate_events = wsp.get_gate_events(sample_id, gate_name=gate, gate_path=gate_path)

    event_labels = get_fluoro_labels(sample, gate_events)
    gate_events = gate_events[event_labels]

    # Calculate the average of the values on the gate events dataframe
    gate_mfi = gate_events.mean()

    # Get the labels from the mfi dataframe
    labels = gate_mfi.index.values.tolist()

    gate = gate.strip()

    if pre == 'ICS':
        # Get the dictionary of indices in order to match the labels to the gates

        # ICS labels
        list1_ics = ['FL1-H IFN-y B525-FITC-H', 'FL1-A IFN-y B525-FITC-A', 'FL3-H TCR GD B690-PC5.5-H',
                     'FL3-A TCR GD B690-PC5.5-A', 'FL4-H TNF-a Y585-PE-H', 'FL4-A TNF-a Y585-PE-A',
                     'FL8-H IL-2 Y763-PC7-H', 'FL8-A IL-2 Y763-PC7-A', 'FL10-H CD8 R712-APCA700-H',
                     'FL10-A CD8 R712-APCA700-A', 'FL11-H CD3 R763-APCA750-H', 'FL11-A CD3 R763-APCA750-A',
                     'FL12-H CD4 V450-PB-H', 'FL12-A CD4 V450-PB-A', 'FL14-H Live Dead V610-H',
                     'FL14-A Live Dead V610-A', 'FL16-H IL-17A V763-H', 'FL16-A IL-17A V763-A']
        list2_ics = ['FL1-H INFg B525-FITC-H', 'FL1-A INFg B525-FITC-A', 'FL8-H CD3 Y763-PC7-H', 'FL8-A CD3 Y763-PC7-A',
                     'FL11-H CD8 R763-APCA750-H', 'FL11-A CD8 R763-APCA750-A', 'FL12-H CD4 V450-PB-H',
                     'FL12-A CD4 V450-PB-A', 'FL14-H Live V610-H', 'FL14-A Live V610-A', 'FL16-H IL17 V763-H',
                     'FL16-A IL17 V763-A']

        list_num = None
        if labels == list1_ics:
            list_num = 1
        elif labels == list2_ics:
            list_num = 2

        indices_dict = labels_dict[list_num]

        idx = indices_dict[gate]

        if gate in indices_dict:  # if the gate is in the dictionary
            mfi = round(gate_mfi[idx], 2)
            if mfi != mfi:  # if there are no events, the mfi dataframe will not have the index found
                mfi_comp[f'{parent} | {gate}'] = 'No Events'
            else:
                mfi_comp[f'{parent} | {gate}'] = mfi

    elif pre == 'Tcell':
        # Get the dictionary of indices in order to match the labels to the gates

        # Tcell labels
        list1_tcell = ['FITC-A', 'BV421-A', 'BV510-A', 'BV605-A', 'BV650-A', 'BV785-A', 'APC-A', 'APC-Alexa 700-A',
                       'APC-Alexa 700-A', 'APC-Cy7-A', 'PE-A', 'PE-Texas Red-A', 'PE-Cy5-5-A', 'PE-Cy7-A']
        # Equal to
        list4_tcell = ['FITC-A CD45RA', 'BV421-A CD57', 'BV510-A CD45', 'BV605-A Live Dead', 'BV650-A LAG-3',
                       'BV711-A CD8intra', 'BV785-A CD95', 'APC-A CD4', 'APC-Alexa 700-A CD8extra',
                       'APC-Alexa 700-A CD8extra', 'APC-Cy7-A CD3', 'PE-A CCR7', 'PE-Texas Red-A CD28', 'PE-Cy7-A CD27',
                       'PerCP-Cy5-5-A PD-1']

        list2_tcell = ['FL1-H CD103 B525-FITC-H', 'FL1-A CD103 B525-FITC-A', 'FL2-H CD39 B610-ECD-H',
                       'FL2-A CD39 B610-ECD-A', 'FL3-H GD B690-PC5.5-H', 'FL3-A GD B690-PC5.5-A',
                       'FL4-H CCR7 Y585-PE-H', 'FL4-A CCR7 Y585-PE-A', 'FL8-H CD69 Y763-PC7-H', 'FL8-A CD69 Y763-PC7-A',
                       'FL9-H 4-1BB R660-APC-H', 'FL9-A 4-1BB R660-APC-A', 'FL10-H CD45RA R712-APCA700-H',
                       'FL10-A CD45RA R712-APCA700-A', 'FL11-H CD3 APC-Cy7 R763-APCA750-H',
                       'FL11-A CD3 APC-Cy7 R763-APCA750-A', 'FL12-H CD4 V450-PB-H', 'FL12-A CD4 V450-PB-A',
                       'FL13-H CD8 V525-KrO-H', 'FL13-A CD8 V525-KrO-A', 'FL14-H LiveDead V610-H',
                       'FL14-A LiveDead V610-A', 'FL15-H LAG3 V660-H', 'FL15-A LAG3 V660-A', 'FL16-H CD95 V763-H',
                       'FL16-A CD95 V763-A']
        # Equal to
        list3_tcell = ['FL1-H B525-FITC-H', 'FL1-A B525-FITC-A', 'FL2-H B610-ECD-H', 'FL2-A B610-ECD-A',
                       'FL3-H B690-PC5.5-H', 'FL3-A B690-PC5.5-A', 'FL4-H Y585-PE-H', 'FL4-A Y585-PE-A',
                       'FL8-H Y763-PC7-H', 'FL8-A Y763-PC7-A', 'FL9-H R660-APC-H', 'FL9-A R660-APC-A',
                       'FL10-H R712-APCA700-H', 'FL10-A R712-APCA700-A', 'FL11-H R763-APCA750-H',
                       'FL11-A R763-APCA750-A', 'FL12-H V450-PB-H', 'FL12-A V450-PB-A', 'FL13-H V525-KrO-H',
                       'FL13-A V525-KrO-A', 'FL14-H V610-H', 'FL14-A V610-A', 'FL15-H V660-H', 'FL15-A V660-A',
                       'FL16-H V763-H', 'FL16-A V763-A']

        list_num = None
        if labels == list1_tcell or labels == list4_tcell:
            list_num = 1
        elif labels == list2_tcell or labels == list3_tcell:
            list_num = 2

        if labels == list1_tcell:
            list1_tcell.insert(5, 'BV711-A')

        indices_dict = labels_dict[list_num]

        gates = {}  # Gates data structure

        # Add gates to the gates data structure given the gate being handled
        if gate in tcm_gate_aliases:
            gates['TCM'] = ['CD45RA', 'CCR7']
        elif gate in tem_gate_aliases:
            gates['TEM'] = ['CD45RA', 'CCR7']
        elif gate in teff_gate_aliases:
            gates['TEFF'] = ['CD45RA', 'CCR7']
        elif gate in precursor_gate_aliases:
            gates['PRECURSOR'] = ['CD45RA', 'CCR7']
        elif len(gate.strip().split(' ')) > 1:
            check_gate_suffix(gate, gates)
        else:
            if gate.endswith('+') or gate.endswith('-'):
                gates[gate] = [gate[:-1]]
            else:
                gates[gate] = [gate]

        extract_gate_mfi(gate_mfi, gates, indices_dict, mfi_comp)

    elif pre == 'Thoming':
        # Get the dictionary of indices in order to match the labels to the gates

        list1_thoming = ['FITC-A CD103', 'PerCP-Cy5-5-A CCR4', 'BV421-A CD4', 'BV510-A CD8', 'BV605-A LIVE',
                         'BV650-A CXCR3', 'BV785-A CD95', 'APC-A CCR9 AF677', 'APC-Alexa 700-A CD45RA',
                         'APC-Alexa 700-A CD45RA', 'APC-Cy7-A CD3', 'PE-A CCR7', 'PE-Cy7-A CCR6']
        list2_thoming = ['FL1-H CD103 B525-FITC-H', 'FL1-A CD103 B525-FITC-A', 'FL2-H CCR4 B610-ECD-H',
                         'FL2-A CCR4 B610-ECD-A', 'FL3-H GD B690-PC5.5-H', 'FL3-A GD B690-PC5.5-A',
                         'FL4-H CCR7 Y585-PE-H', 'FL4-A CCR7 Y585-PE-A', 'FL8-H CCR6 Y763-PC7-H',
                         'FL8-A CCR6 Y763-PC7-A', 'FL9-H CCR9 R660-APC-H', 'FL9-A CCR9 R660-APC-A',
                         'FL10-H CD45RA R712-APCA700-H', 'FL10-A CD45RA R712-APCA700-A', 'FL11-H CD3 R763-APCA750-H',
                         'FL11-A CD3 R763-APCA750-A', 'FL12-H CD4 V450-PB-H', 'FL12-A CD4 V450-PB-A',
                         'FL13-H CD8 V525-KrO-H', 'FL13-A CD8 V525-KrO-A', 'FL14-H Live dead V610-H',
                         'FL14-A Live dead V610-A', 'FL15-H CxCR3 V660-H', 'FL15-A CxCR3 V660-A', 'FL16-H CD95 V763-H',
                         'FL16-A CD95 V763-A']
        list3_thoming = ['FITC-A CD103', 'BV421-A CD4', 'BV510-A CD8', 'BV605-A Live Dead', 'BV650-A CXCR3',
                         'BV711-A CD8 intra', 'BV785-A CD95', 'APC-A CCR9', 'APC-Alexa 700-A CD45RA',
                         'APC-Alexa 700-A CD45RA', 'APC-Cy7-A CD3', 'PE-A CCR7', 'PE-Texas Red-A CCR4', 'PE-Cy7-A CCR6',
                         'PerCP-Cy5-5-A TCRgd']
        list4_thoming = ['FL1-H CD103 B525-FITC-H', 'FL1-A CD103 B525-FITC-A', 'FL2-H CCR4 B610-ECD-H',
                         'FL2-A CCR4 B610-ECD-A', 'FL4-H CCR7 Y585-PE-H', 'FL4-A CCR7 Y585-PE-A',
                         'FL7-H TCRgd Y710-PC5.5-H', 'FL7-A TCRgd Y710-PC5.5-A', 'FL8-H CCR6 Y763-PC7-H',
                         'FL8-A CCR6 Y763-PC7-A', 'FL9-H CCR9 R660-APC-H', 'FL9-A CCR9 R660-APC-A',
                         'FL10-H CD45RA R712-APCA700-H', 'FL10-A CD45RA R712-APCA700-A', 'FL11-H CD3 R763-APCA750-H',
                         'FL11-A CD3 R763-APCA750-A', 'FL12-H CD4 V450-PB-H', 'FL12-A CD4 V450-PB-A',
                         'FL13-H CD8 V525-KrO-H', 'FL13-A CD8 V525-KrO-A', 'FL14-H Live dead V610-H',
                         'FL14-A Live dead V610-A', 'FL15-H CXCR3 V660-H', 'FL15-A CXCR3 V660-A', 'FL16-H CD95 V763-H',
                         'FL16-A CD95 V763-A']

        list_num = None
        if labels == list1_thoming:
            list_num = 1
        elif labels == list2_thoming:
            list_num = 2
        elif labels == list3_thoming:
            list_num = 3
        elif labels == list4_thoming:
            list_num = 4

        indices_dict = labels_dict[list_num]

        gates = {}  # Gates data structure

        # Add gates to the gates data structure given the gate being handled
        if th1_gate(gate, gate_path, parent):
            gates['Th1'] = ['CXCR3', 'CCR6']
        elif th2_gate(gate, gate_path, parent):
            gates['Th2'] = ['CCR4', 'CCR6']
        elif th17_gate(gate, gate_path, parent):
            gates['Th17'] = ['CCR4', 'CCR6']
        elif th1star_gate(gate, gate_path, parent):
            gates['Th1*'] = ['CCR4', 'CXCR3', 'CCR6']
        elif tscm_gate(gate, gate_path, parent):
            gates['Tscm'] = ['CD45RA', 'CCR7', 'CD95']
        elif len(gate.strip().split(' ')) > 1:
            check_gate_suffix(gate, gates)
        else:
            if gate.endswith('+') or gate.endswith('-'):
                gates[gate] = [gate[:-1]]
            else:
                gates[gate] = [gate]

        extract_gate_mfi(gate_mfi, gates, indices_dict, mfi_comp)


def extract_gate_mfi(gate_mfi, gates, indices_dict, mfi_comp):
    """
        Get the MFI for the gate or gates from the mfi dataframe and store it

        :param gate_mfi: the MFI dataframe
        :param gates: the gates to extract
        :param indices_dict: the dictionary for the indices where all gates are located in the mfi dataframe
        :param mfi_comp: the MFI comparison data structure
        :return: None
    """

    for label, gate_list in gates.items():
        if type(gate_list) == list:
            mfi_comp[label] = {}
            for item in gate_list:
                idx = indices_dict[item]
                mfi = float(round(gate_mfi[idx], 2))
                if mfi != mfi:
                    mfi_comp[label][item] = 'No Events'
                else:
                    mfi_comp[label][item] = mfi
        else:
            idx = indices_dict[gate_list]
            mfi = round(gate_mfi[idx], 2)
            if mfi != mfi:
                mfi_comp[label] = 'No Events'
            else:
                mfi_comp[label] = mfi


def check_gate_suffix(gate, gates):
    """
        Checks the gate suffix in order to match the indices dictionary
        :param gate: the gate being handled
        :param gates: the gates data structure
        :return: None
    """

    gate = gate.strip()
    gates[gate] = []
    for item in gate.split(' '):
        if item != '+' and item != '-':
            if item[-1] == '+' or item[-1] == '-':
                gates[gate].append(item[:-1])
            else:
                gates[gate].append(item)


def th1star_gate(gate, gate_path, parent):
    """
        Checks if the gate being handled is the Th1* gate
        :param gate: the gate being handled
        :param gate_path: the gate path
        :param parent: the parent of the gate
        :return: true if it is, false otherwise
    """

    return 'th1+' in str(gate).lower() or 'th1*' in str(gate).lower() or (
            (gate == 'CCR4+ CXCR3+' and parent == 'CCR6+') and 'CD4+' in gate_path)


def th17_gate(gate, gate_path, parent):
    """
        Checks if the gate being handled is the Th17 gate
        :param gate: the gate being handled
        :param gate_path: the gate path
        :param parent: the parent of the gate
        :return: true if it is, false otherwise
    """

    return 'th17' in str(gate).lower() or ((gate == 'CCR4+' and parent == 'CCR6+') and 'CD4+' in gate_path)


def th2_gate(gate, gate_path, parent):
    """
        Checks if the gate being handled is the Th2 gate
        :param gate: the gate being handled
        :param gate_path: the gate path
        :param parent: the parent of the gate
        :return: true if it is, false otherwise
    """

    return 'th2' in str(gate).lower() or ((gate == 'CCR4+' and parent == 'CCR6-') and 'CD4+' in gate_path)


def th1_gate(gate, gate_path, parent):
    """
        Checks if the gate being handled is the Th1 gate
        :param gate: the gate being handled
        :param gate_str: the gate converted to lowercase
        :param gate_path: the gate path
        :param parent: the parent of the gate
        :return: true if it is, false otherwise
    """

    th1 = 'th1' in str(gate).lower() and 'th17' not in str(gate).lower() and 'th1+' not in str(
        gate).lower() and 'th1*' not in str(gate).lower()
    return th1 or ((gate == 'CXCR3+' and parent == 'CCR6-') and 'CD4+' in gate_path)


def tscm_gate(gate, gate_path, parent):
    """
        Checks if the gate being handled is the Tscm gate
        :param gate: the gate being handled
        :param gate_str: the gate being handled in converted to lowercase
        :param gate_path: the gate path
        :param parent: the parent of the gate
        :return: true if it is, false otherwise
    """

    return 'tscm' in str(gate).lower() or ((gate == 'CD45RA+ CCR7+' and parent == 'CD95+') and 'CD8+' in gate_path)


def thoming_gate_of_interest(gate, gate_path):
    """
        Checks if the gate is one of the Thoming gates to be analyzed
        :param gate: the gate being handled
        :param gate_path: the gate path
        :return: true if it is, false otherwise
    """

    gate_str = str(gate).lower()
    parent = gate_path[-1]
    th1_gate = 'th1' in gate_str and 'th17' not in gate_str and 'th1+' not in gate_str and 'th1*' not in gate_str
    tscm_gate = 'tscm' in gate_str or ((gate == 'CD45RA+ CCR7+' and parent == 'CD95+') and 'CD8+' in gate_path)
    return tscm_gate or \
        (((th1_gate or 'th2' in gate_str or 'th1*' in gate_str or 'th1+' in gate_str or 'th17' in gate_str) or
          (gate == 'CCR4+ CXCR3+' and parent == 'CCR6+') or  # Th1*
          (gate == 'CCR4+' and parent == 'CCR6+') or  # Th17
          (gate == 'CCR4+' and parent == 'CCR6-') or  # Th2
          (gate == 'CXCR3+' and parent == 'CCR6-')) and 'CD4+' in gate_path)  # Th1


def write_sample_id(sample_id):
    with open('Samples/sample_ids.txt', 'a') as f:
        f.write(f'{sample_id[:-4]}\n')


def analyze(wsp_file, analysis, pre):
    """
        Analyze the workspace and its corresponding sample files given the type of
        analysis(TIL, PBMC) and the prefix (ICS, Tcell, Thoming)

        :param wsp_file: The workspace file path
        :param analysis: The type of analysis
        :param pre: The prefix of the analysis

        :return: None
    """

    print(f'\nWORKSPACE {wsp_file}\n')

    # Create a Workspace with the path to our WSP file and FCS files
    wsp = fk.Workspace(wsp_file, fcs_samples=os.path.dirname(wsp_file), ignore_missing_files=True)

    # Loop through all sample groups
    sample_group = 'All Samples'

    # Analyze samples in order to fetch analysis results
    wsp.analyze_samples(sample_group)

    # Get sample file names
    sample_list = wsp.get_sample_ids(group_name=sample_group)

    # loop through samples
    for sample_id in sample_list:
        print("##################################################")
        print(f'\nSAMPLE {sample_id}\n')

        write_sample_id(sample_id)

        # Get the sample from its ID
        sample = wsp.get_sample(sample_id)

        # Get the sample date
        date = sample.acquisition_date

        # Get gate hierarchy
        hierarchy = wsp.get_gate_hierarchy(sample_id)
        print(hierarchy)

        # Get the sample results dataframe
        sample_results = wsp.get_gating_results(sample_id)

        # Data structures to hold our data
        mfi_comp = {}
        gate_pct = {}

        # gate aliases checking for Gamma Delta + gate
        gd_gate_name = ''
        for alias in gd_gate_aliases:
            if alias in hierarchy:
                gd_gate_name = alias
                break

        if pre == 'ICS':
            # gate aliases checking for IFN-g gate
            ifng_gate_name = ''
            for alias in ifng_gate_aliases:
                if alias in hierarchy:
                    ifng_gate_name = alias
                    break

            ics_gates = [f'CD4+/{ifng_gate_name}', f'CD8+/{ifng_gate_name}', f'{gd_gate_name}/{ifng_gate_name}']

            labels_dict_ics = {
                1: {ifng_gate_name: 1, 'CD3+': 11, 'CD4+': 13, 'CD8+': 9, gd_gate_name: 3},
                2: {ifng_gate_name: 1, 'CD3+': 3, 'CD4+': 7, 'CD8+': 5}
            }

            # The gates we want to extract data from
            for item in ics_gates:
                parent, gate = item.split('/')
                # The gate structure
                for gate_id, gate_path in wsp.get_gate_ids(sample_id):
                    # if the gate id is a gate of interest, the last item in the gate path is the parent
                    # and the gate is a descendant of the CD3+ gate
                    if gate_id == gate and gate_path[-1] == parent:
                        # extract_data function
                        extract_data(gate, gate_path, gate_pct, mfi_comp, parent, sample, sample_id,
                                     sample_results, wsp, labels_dict_ics)
        elif pre == 'Tcell':
            # gate aliases checking for LAG-3 gate
            lag3_gate_name = ''
            for alias in lag3_gate_aliases:
                if alias in hierarchy:
                    lag3_gate_name = alias
                    break

            # gate aliases checking for TCM gate
            tcm_gate_name = ''
            for alias in tcm_gate_aliases:
                if alias in hierarchy:
                    tcm_gate_name = alias
                    break

            # gate aliases checking for TEM gate
            tem_gate_name = ''
            for alias in tem_gate_aliases:
                if alias in hierarchy:
                    tem_gate_name = alias
                    break

            # gate aliases checking for Teff gate
            teff_gate_name = ''
            for alias in teff_gate_aliases:
                if alias in hierarchy:
                    teff_gate_name = alias
                    break

            # gate aliases checking for Precursor gate
            precursor_gate_name = ''
            for alias in precursor_gate_aliases:
                if alias in hierarchy:
                    precursor_gate_name = alias
                    break

            if lag3_gate_name.endswith('+'):
                lag3_gate_name_plus = lag3_gate_name[:-1]
            else:
                lag3_gate_name_plus = lag3_gate_name

            labels_dict_tcell = {
                1: {'CD45RA': 0, 'CD57': 1, 'CD45': 2, precursor_gate_name: 3, lag3_gate_name_plus.strip(): 4,
                    'CD95': 6, 'CD4': 7, 'CD8': 8,
                    'CD3': 10, 'CCR7': 11, 'CD28': 12, 'CD27': 13, 'PD-1': 14},
                2: {'CD103': 1, 'CD39': 3, 'GD': 5, 'CCR7': 7, 'CD69': 9, '4-1BB': 11, 'CD45RA': 13, 'CD3': 15,
                    'CD4': 17, 'CD8': 19, precursor_gate_name: 21, lag3_gate_name_plus.strip(): 23, 'CD95': 25},
            }

            if analysis == 'PBMCs':
                # Check if the CD4+ gate (alone) is in the gating strategy.
                # If it is not, it is joined with the CD8+ gate
                cd4_in = False
                for gate, _ in wsp.get_gate_ids(sample_id):
                    if gate == 'CD4+':
                        cd4_in = True
                        break

                if not cd4_in:
                    tcell_pbmc_gates = [f'CD4+ CD8+/{precursor_gate_name}', f'CD4+ CD8+/{tcm_gate_name}',
                                        f'CD4+ CD8+/{tem_gate_name}', f'CD4+ CD8+/{teff_gate_name}',
                                        'CD4+ CD8+/CD95+', f'CD4+ CD8+/{lag3_gate_name}', 'CD4+ CD8+/CD103+',
                                        f'GD+/{precursor_gate_name}', f'GD+/{tcm_gate_name}', f'GD+/{tem_gate_name}',
                                        f'GD+/{teff_gate_name}', 'GD+/CD95+', f'GD+/{lag3_gate_name}', 'GD+/CD103+']
                else:
                    tcell_pbmc_gates = [f'CD4+/{tcm_gate_name}',
                                        f'CD4+/{tem_gate_name}', f'CD4+/{teff_gate_name}',
                                        'CD4+/CD95+', f'CD4+/{lag3_gate_name}', 'CD4+/CD103+',
                                        f'CD8+/{tcm_gate_name}',
                                        f'CD8+/{tem_gate_name}', f'CD8+/{teff_gate_name}',
                                        'CD8+/CD95+', f'CD8+/{lag3_gate_name}', 'CD8+/CD103+', f'GD+/{tcm_gate_name}',
                                        f'GD+/{tem_gate_name}',
                                        f'GD+/{teff_gate_name}', 'GD+/CD95+', f'GD+/{lag3_gate_name}', 'GD+/CD103+']

                # The gates we want to extract data from
                for item in tcell_pbmc_gates:
                    parent, gate = item.split('/')
                    # The gate structure
                    for gate_id, gate_path in wsp.get_gate_ids(sample_id):
                        # if the gate id is a gate of interest, the last item in the gate path is the parent
                        # and the gate is a descendant of the CD3+ gate
                        if gate in gate_id and gate_path[-1] == parent and 'CD3+' in gate_path \
                                or str(gate) == 'CD39- CD69-':
                            extract_data(gate, gate_path, gate_pct, mfi_comp, parent, sample,
                                         sample_id, sample_results, wsp, labels_dict_tcell)

            elif analysis == 'TILs':
                # The gate structure
                for gate, gate_path in wsp.get_gate_ids(sample_id):
                    # Use everything except (CD4+ CD8+), all gates that have 2 minus (-) signs
                    # and all gates with no child gates
                    if gate != 'CD4+ CD8+' and 'CD4+ CD8+' not in gate_path and str(gate).count('-') != 2 \
                            and not wsp.get_child_gate_ids(sample_id, gate, gate_path) or str(gate) == 'CD39- CD69-':
                        parent = gate_path[-1]
                        extract_data(gate, gate_path, gate_pct, mfi_comp, parent, sample,
                                     sample_id, sample_results, wsp, labels_dict_tcell)

        elif pre == 'Thoming':

            # gate aliases checking for Precursor gate
            precursor_gate_name = ''
            for alias in precursor_gate_aliases:
                if alias in hierarchy:
                    precursor_gate_name = alias
                    break

            labels_dict_thoming = {
                1: {'CD103': 0, 'CCR4': 1, 'CD4': 2, 'CD8': 3, precursor_gate_name: 4, 'CXCR3': 5, 'CD95': 6,
                    'CCR9': 7, 'CD45RA': 8, 'CD3': 10, 'CCR7': 11, 'CCR6': 12},
                2: {'CD103': 1, 'CCR4': 3, gd_gate_name: 5, 'CCR7': 7, 'CCR6': 9, 'CCR9': 11, 'CD45RA': 13, 'CD3': 15,
                    'CD4': 17, 'CD8': 19, precursor_gate_name: 21, 'CXCR3': 23, 'CD95': 25},
                3: {'CD103': 0, 'CD4': 1, 'CD8': 2, precursor_gate_name: 3, 'CXCR3': 4, 'CD95': 6, 'CCR9': 7,
                    'CD45RA': 8, 'CD3': 10, 'CCR7': 11, 'CCR4': 12, 'CCR6': 13, gd_gate_name: 14},
                4: {'CD103': 1, 'CCR4': 3, 'CCR7': 5, gd_gate_name: 7, 'CCR6': 9, 'CCR9': 11, 'CD45RA': 13,
                    'CD3': 15, 'CD4': 17, 'CD8': 19, precursor_gate_name: 21, 'CXCR3': 23, 'CD95': 25}
            }

            # The gate structure
            for gate, gate_path in wsp.get_gate_ids(sample_id):
                # Check for the gates of interest
                if thoming_gate_of_interest(gate, gate_path):
                    extract_data(gate.strip(), gate_path, gate_pct, mfi_comp, gate_path[-1], sample,
                                 sample_id, sample_results, wsp, labels_dict_thoming)

        print(mfi_comp)
        print(gate_pct)

        # Cut off the file type extension
        sample_id_path = sample_id[: -4]

        # Create the directory structure
        if not os.path.exists(f'Samples/{pre}_{analysis}/{sample_id_path}_{date}'):
            os.makedirs(f'Samples/{pre}_{analysis}/{sample_id_path}_{date}')

        # Build directed graph to extract gate hierarchy
        dot = graphviz.Digraph(f'{sample_id_path}_gate_hierarchy', strict=True)
        for gate, gate_path in wsp.get_gate_ids(sample_id):
            if ':' in gate:
                gate = gate.replace(':', '')
            dot.node(gate)
            for i in range(len(gate_path)):
                if i == len(gate_path) - 1:
                    dot.edge(gate_path[i], gate)
                else:
                    dot.edge(gate_path[i], gate_path[i + 1])
        dot.format = 'png'
        dot.render(directory=f'Samples/{pre}_{analysis}/{sample_id_path}_{date}') \
            .replace('\\', '/')

        # if the mfi dataframe is not empty (i.e. we found gates of interest)
        if mfi_comp:
            with open(f'{os.getcwd()}/Samples/{pre}_{analysis}/{sample_id_path}_{date}/{sample_id_path}_mfi.csv', 'w',
                      newline='\n') as f:
                writer = csv.DictWriter(f, fieldnames=mfi_comp.keys())
                writer.writeheader()
                writer.writerow(mfi_comp)

            with open(f'{os.getcwd()}/Samples/{pre}_{analysis}/global_mfi.csv', 'a', newline='\n') as f_object:
                writer_object = csv.DictWriter(f_object, fieldnames=mfi_comp.keys())
                writer_object.writeheader()
                writer_object.writerow(mfi_comp)

        # if the gate percentages dictionary is not empty (i.e. we found gates of interest)
        if gate_pct:
            with open(f'Samples/{pre}_{analysis}/{sample_id_path}_{date}/{sample_id_path}_gate_percentages.csv', 'w',
                      newline='\n') as f:
                writer = csv.DictWriter(f, fieldnames=gate_pct.keys())
                writer.writeheader()
                writer.writerow(gate_pct)

            with open(f'{os.getcwd()}/Samples/{pre}_{analysis}/global_gate_pct.csv', 'a', newline='\n') as f_object:
                writer_object = csv.DictWriter(f_object, fieldnames=gate_pct.keys())
                writer_object.writeheader()
                writer_object.writerow(gate_pct)
        print("##################################################")


# ---------------------------------------------------------------------------------
# Execution starts here

# base directory where the files are stored
base_dir = "Data/DATA_Raw_files/FLOW"

# All gate aliases found
gd_gate_aliases = ['Gamma delta + ', 'GD+', 'TCRgd+', 'TCR gd+', 'gd+']
ifng_gate_aliases = ['INFg+', 'IFN-g+']
lag3_gate_aliases = ['LAG3+', 'LAG3', 'LAG-3 +', 'Lag3+']
tcm_gate_aliases = ['TCM CD45RA- CCR7+', 'TCM']
teff_gate_aliases = ['Teff CD45RA+ CCR7-', 'Teff', 'TEFF']
tem_gate_aliases = ['TEM CD45RA- CCR7-', 'TEM']
precursor_gate_aliases = ['precursors', 'precursor', 'Precursor', 'PRECURSOR', 'Naive CD45RA+ CCR7+', 'Live Dead',
                          'Live-Dead']

if __name__ == '__main__':

    if os.path.exists(f'{os.getcwd()}/Samples'):
        shutil.rmtree(f'{os.getcwd()}/Samples', ignore_errors=True)

    os.mkdir('Samples')

    wsp_files = get_wsp_files(base_dir)

    # Loop through workspace files
    for wsp_file in wsp_files:

        analysis = ''
        if 'PBMC' in wsp_file:
            analysis = 'PBMCs'
        elif 'TIL' in wsp_file:
            analysis = 'TILs'

        pre = ''
        if 'ICS' in wsp_file:
            pre = 'ICS'
        elif 'Tcell' in wsp_file:
            pre = 'Tcell'
        elif 'Thoming' in wsp_file:
            pre = 'Thoming'

        analyze(wsp_file, analysis, pre)
