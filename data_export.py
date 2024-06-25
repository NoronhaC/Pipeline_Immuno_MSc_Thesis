import mysql.connector
import tkinter as tk
from tkinter import messagebox
import ttkbootstrap as ttk


# Functions

def query_database(query):
    res = None
    cursor.execute(query)
    if query.lower().startswith('select') or query.lower().startswith('show'):
        res = cursor.fetchall()
    return res


def build_form(master_window, table):
    cols = table_cols[table]
    print(cols)
    for col in cols:
        print(f'{col} : {form_items[col]}')
    pass


def build_menu(table):
    if len(table) == 0:
        messagebox.showerror('Table Undefined', 'Please choose a table to insert data!')
    else:
        new_window = ttk.Toplevel(root)
        menu_frame = tk.Frame(master=new_window, background='red')
        menu_frame.pack(fill="both", expand=True, padx=20, pady=20)
        build_form(new_window, table, table_cols[table])


def change_label(menu, txt):
    menu['text'] = txt


# ---------------------------------------------------------------------------------
# Execution starts here

tables = ['elisa', 'output', 'patient', 'phenotype_cytometry', 'sample', 'tumor_biopsy',
        'wes_whole_exosome_sequencing', 'whole_blood']

table_cols = {
    'elisa': ('Patient_ID', 'Sample_ID', 'Protein_Cytokine_Interest', 'ELISA_Date', 'ELISA_Plate_ID'),
    'output': ('Output_ID', 'Patient_ID', 'Sample_ID', 'ELISA_CytokineProtein_Conc_Peptide', 'Flow_Panel_ID',
               'CD_Percentage', 'CD_MFI', 'CytokineProtein_Percentage', 'CytokineProtein_MFI'),
    'patient': ('Patient_ID', 'Age', 'Sex'),
    'phenotype_cytometry': ('Flow_Panel_ID', 'Patitent_ID', 'Sample_ID'),
    'sample': ('Sample_ID', 'Patient_ID', 'Tumor_Biopsy_ID', 'Whole_Blood_ID', 'Tissue_Type', 'Sample_Date'),
    'tumor_biopsy': ('Tumor_Biopsy_ID', 'WES_Whole_Exosome_Sequencing', 'Biopsy_Location', 'Tissue_Type',
                     'Pathology_Diagnosis'),
    'wes_whole_exosome_sequencing': ('Artificial_Peptides_from', 'Biopsy_Location', 'Tissue_Type', 'Pathology_Diagnosis'),
    'whole_blood': ('Whole_Blood_ID', 'Heparin', 'EDTA', 'Serum', 'LRS_chamber', 'Buffy_Coat')
}

flow_panel_id = ['ICS', 'T-cell panel', 'Thoming']

pathology_diagnosis = ['Intrahepatic Cholangiocarcinoma', 'Metastatic Colon Adenocarcinoma',
                       'Colon Adenocarcinoma', 'Colon Neuroendocrine Tumor', 'Rectum Adenocarcinoma',
                       'Pancreas Neuroendocrine Tumor', 'Pancreas Adenocarcinoma', 'Retroperitoneum Liposarcoma',
                       'Metastatic Melanoma possibly, uveal melanoma',
                       'Glioblastoma', 'Ductal Adenocarcinoma',
                       'Ampullary Adenomyomatous Hyperplasia',
                       'Stomach Gastrointestinal Stromal Tumour',
                       'Billiary Tract Adenocarcinoma']


def get_samples():
    with open('Samples/sample_ids.txt', 'r') as f:
        lines = f.readlines()
        lines = set(list(map(lambda x: x[:-1], lines)))
        return lines


def get_peptides():
    with open('Patients/peptides.txt', 'r') as f:
        lines = f.readlines()
        lines = set(list(map(lambda x: x[:-1], lines)))
        return lines


form_items = {
    'Patient_ID': {'label': 'Select Patient:', 'data_type': 'dropdown', 'dropdown_list': [1, 2, 3, 4, 5]},
    'Sample_ID': {'label': 'Select Sample:', 'data_type': 'dropdown', 'dropdown_list': get_samples()},
    'Protein_Cytokine_Interest': {'label': 'Select Cytokine:', 'data_type': 'str'},
    'ELISA_Date': {'label': 'Select Date:', 'data_type': 'date'},
    'ELISA_Plate_ID': {'label': 'Select Plate:', 'data_type': 'int'},
    'Output_ID': {'label': '', 'data_type': 'inc_int'},
    'ELISA_CytokineProtein_Conc_Peptide': {'label': 'Select Peptide', 'data_type': 'dropdown', 'dropdown_list': get_peptides()},
    'Flow_Panel_ID': {'label': 'Select Flow Panel ID:', 'data_type': 'dropdown', 'dropdown_list': flow_panel_id},
    'Pathology_Diagnosis': {'label': 'Select Pathology Diagnosis:', 'data_type': 'dropdown', 'dropdown_list': pathology_diagnosis},
}

root = ttk.Window(themename='superhero')
root.title("User Input")
root.resizable(False, False)

main_frame = ttk.Frame(master=root)
main_frame.pack(fill="both", expand=True, padx=20, pady=20)

table_var_label = ttk.StringVar()

table_label = ttk.Label(main_frame, textvariable=table_var_label)
table_var_label.set('Table Name')
table_label.grid(row=0, column=1, padx=5)

# Flow Panel menu
clicked_table_var = ttk.StringVar()

dropdown_table = ttk.Menubutton(main_frame, bootstyle='primary', text='Pick the Table')
dropdown_table.grid(row=0, column=2)

inside_dropdown_table = ttk.Menu(dropdown_table)
for x in tables:
    inside_dropdown_table.add_radiobutton(label=x, variable=clicked_table_var,
                                          command=lambda x=x: change_label(dropdown_table, x))

dropdown_table['menu'] = inside_dropdown_table

next_button = ttk.Button(root, text="Next", command=lambda: build_menu(clicked_table_var.get()),
                         bootstyle='success, outline')
next_button.pack(pady=10)

root.mainloop()

# # Pathology Diagnosis menu
# clicked_pathology_diagnosis = ttk.StringVar()
#
# dropdown_path_diag = ttk.Menubutton(root, bootstyle='primary', text='Pick the Pathology Diagnosis')
# dropdown_path_diag.pack(pady=10)
#
# inside_dropdown_path_diag = ttk.Menu(dropdown_path_diag)
# for x in pathology_diagnosis:
#     inside_dropdown_path_diag.add_radiobutton(label=x, variable=clicked_pathology_diagnosis,
#     command=lambda x=x: change_label(dropdown_path_diag, x))
#
# dropdown_path_diag['menu'] = inside_dropdown_path_diag
#
# # Button for closing
# exit_button = ttk.Button(root, text="Save & Exit", command=root.destroy, bootstyle='success, outline')
# exit_button.pack(pady=10, side=BOTTOM)

conn = mysql.connector.connect(
    host='localhost',
    user='root',
    password='root',
    database='biomedicaldb'
)

cursor = conn.cursor()

# tables = [table[0] for table in query_database('SHOW TABLES')]
# print(tables)
#
# for table in tables:
#     entries = query_database(f'SELECT * FROM biomedicaldb.{table}')
#     print(f'{table} - {cursor.column_names}')
