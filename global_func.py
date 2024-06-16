import os


# Functions

def read_csv(csv_file):
    """
        Read the data from a csv and transform it into a dictionary
    :param csv_file: the csv file path
    :return: the data from the csv file path, in a dictionary
    """
    csv_data = []
    with open(csv_file, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            idx = lines.index(line)
            if idx % 2 == 1:
                if 'mfi' in csv_file and 'ICS' not in csv_file:
                    values_list = [x for x in line.split('\"') if x.strip() != '' and x != ',']
                    copy_values_list = values_list[:]
                    for item in copy_values_list:
                        if item.startswith(',') or item.endswith(','):
                            item_idx = values_list.index(item)
                            values_list.remove(item)
                            aux_list = [x for x in item.split(',') if x != '']
                            for i in range(len(aux_list)):
                                values_list.insert(item_idx + i, aux_list[i])
                    csv_data.append(dict(zip(lines[idx - 1].split(','), values_list)))
                else:
                    csv_data.append(dict(zip(lines[idx - 1].split(','), line.split(','))))
    return csv_data


def get_data(prefix, group, analysis, data_type):
    """
        Get the data given the conditions
    :param prefix: FLOW or ELISA
    :param group: the group in which the prefix divides
    :param analysis: PBMCs, TILs or WBA
    :param data_type: gate mfi, gate percentages
    :return: the data, as a list
    """
    if prefix.upper() == 'ELISA':
        return data[prefix][group][analysis]
    else:
        if group == '':
            result = []
            for aux_group in ['ICS', 'Tcell', 'Thoming']:
                result.extend(data[prefix][aux_group][analysis][data_type])
            return result
        return data[prefix][group][analysis][data_type]


def get_gate_data(gate, prefix, group, analysis, data_type):
    """
        Get data from a specific gate
    :param gate: the gate
    :param prefix: FLOW
    :param group: ICS, Tcell or Thoming
    :param analysis: PBMCs or TILs
    :param data_type: gate percentages or gate mfi
    :return: the data, as a list
    """
    res_data = []
    gate_data = get_data(prefix, group, analysis, data_type)
    for dict_item in gate_data:
        for key, value in dict_item.items():
            if gate in key.lower():
                res_data.append({key: value})
    return res_data


# ---------------------------------------------------------------------------------
# Execution starts here

# Import Data

# Build data structure
data = {'ELISA': {'Concentrations': {'TIL': [], 'WBA': []}, 'OD': {'TIL': [], 'WBA': []}, 'Reactions': {'TIL': [], 'WBA': []}},
        'FLOW': {'ICS': {'PBMCs': {'gate_mfi': {}, 'gate_pct': {}}, 'TILs': {'gate_mfi': {}, 'gate_pct': {}}},
                 'Tcell': {'PBMCs': {'gate_mfi': {}, 'gate_pct': {}}, 'TILs': {'gate_mfi': {}, 'gate_pct': {}}},
                 'Thoming': {'PBMCs': {'gate_mfi': {}, 'gate_pct': {}}, 'TILs': {'gate_mfi': {}, 'gate_pct': {}}}
                 }
        }

# Import ELISA data
data['ELISA']['Concentrations']['TIL'] = read_csv('Patients/TIL_global_concentrations.csv')
data['ELISA']['Concentrations']['WBA'] = read_csv('Patients/WBA_global_concentrations.csv')
data['ELISA']['OD']['TIL'] = read_csv('Patients/TIL_global_ODs.csv')
data['ELISA']['OD']['WBA'] = read_csv('Patients/WBA_global_ODs.csv')
data['ELISA']['Reactions']['TIL'] = read_csv('Patients/TIL_global_peptide_reactions.csv')
data['ELISA']['Reactions']['WBA'] = read_csv('Patients/WBA_global_peptide_reactions.csv')

# Import FLOW data
for directory in os.listdir('Samples'):
    for item in os.listdir(os.path.join('Samples', directory)):
        path = os.path.join('Samples', directory, item)
        if os.path.isfile(path):
            if 'ICS_PBMCs' in path:
                if 'mfi' in path:
                    data['FLOW']['ICS']['PBMCs']['gate_mfi'] = read_csv(path)
                elif 'gate_pct' in path:
                    data['FLOW']['ICS']['PBMCs']['gate_pct'] = read_csv(path)
            elif 'ICS_TILs' in path:
                if 'mfi' in path:
                    data['FLOW']['ICS']['TILs']['gate_mfi'] = read_csv(path)
                elif 'gate_pct' in path:
                    data['FLOW']['ICS']['TILs']['gate_pct'] = read_csv(path)
            elif 'Tcell_PBMCs' in path:
                if 'mfi' in path:
                    data['FLOW']['Tcell']['PBMCs']['gate_mfi'] = read_csv(path)
                elif 'gate_pct' in path:
                    data['FLOW']['Tcell']['PBMCs']['gate_pct'] = read_csv(path)
            elif 'Tcell_TILs' in path:
                if 'mfi' in path:
                    data['FLOW']['Tcell']['TILs']['gate_mfi'] = read_csv(path)
                elif 'gate_pct' in path:
                    data['FLOW']['Tcell']['TILs']['gate_pct'] = read_csv(path)
            elif 'Thoming_PBMCs' in path:
                if 'mfi' in path:
                    data['FLOW']['Thoming']['PBMCs']['gate_mfi'] = read_csv(path)
                elif 'gate_pct' in path:
                    data['FLOW']['Thoming']['PBMCs']['gate_pct'] = read_csv(path)
            elif 'Thoming_TILs' in path:
                if 'mfi' in path:
                    data['FLOW']['Thoming']['TILs']['gate_mfi'] = read_csv(path)
                elif 'gate_pct' in path:
                    data['FLOW']['Thoming']['TILs']['gate_pct'] = read_csv(path)
