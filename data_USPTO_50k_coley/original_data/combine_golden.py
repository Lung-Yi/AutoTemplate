import pandas as pd

input_file = "data_processed.csv"
input_df = pd.read_csv(input_file)
output_lines = []

for i in range(len(input_df)):
    rxn_smiles = input_df.iloc[i]["rxn_smiles"]
    reaction_id = input_df.iloc[i]["id"]
    line = rxn_smiles + "\t" + reaction_id + "_" + str(i) + "\n"
    output_lines.append(line)

with open("../USPTO_50k_coley.txt", "w") as g:
    g.writelines(output_lines)