"""
Get the reaction ID of each .xlsx file, and output a .txt file in the reaction directory
"""
import os
import pandas as pd
from tqdm import tqdm

input_data_path = "../data_reaxys"
query_template = "{\"fileName\":\"#####\",\"version\":\"1.0\",\"content\":{\"id\":\"root\",\"facts\":[{\"id\":\"Reaxys487\",\"fields\":[{\"value\":\"$$$$$\",\"boundOperator\":\"op_num_equal\",\"id\":\"RX.ID\",\"displayName\":\"Reaction ID\"}],\"fieldsLogicOperator\":\"AND\",\"exist\":false,\"bio\":false}],\"exist\":false,\"bio\":false}}"

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
            for i in range(len(reaction_id_set) // 5000 + 1):
                reaction_ids = ";".join(reaction_id_set[i*5000: (i+1)*5000])
                file_name = f"{reaction_type}_query_{i}.json"
                query_content = query_template.replace("#####", file_name)
                query_content = query_content.replace("$$$$$", reaction_ids)

                with open(os.path.join(reaction_dtr, file_name), "w") as f:
                    f.write(query_content)