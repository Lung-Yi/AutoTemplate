"""
Get the reaction ID of each .xlsx file, and output a .txt file in the reaction directory
"""
import os
import pandas as pd
from tqdm import tqdm

input_data_path = "../data_reaxys"

for reaction_type in tqdm(os.listdir(input_data_path)):
    reaction_dtr = os.path.join(input_data_path, reaction_type)
    reaction_id_set = set()
    for root, dirs, files in os.walk(reaction_dtr, topdown=False):
        for name in files:
            if name.endswith(".xlsx"):
                xlsx_file_path = os.path.join(root, name)
                print(xlsx_file_path)
                excel_data = pd.read_excel(xlsx_file_path, engine='openpyxl')
                excel_data = excel_data[:-3]
                reaction_id_set.update(set(excel_data["Reaction ID"]))
            else:
                pass
        else:
            reaction_id_set = sorted(list(reaction_id_set))
            reaction_id_set = [rxn_id + "\n" for rxn_id in reaction_id_set]
            with open(os.path.join(reaction_dtr, "{}_reaction_id.txt".format(reaction_type)), "w") as f:
                f.writelines(reaction_id_set)
