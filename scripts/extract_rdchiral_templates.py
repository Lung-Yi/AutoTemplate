from rdchiral.template_extractor import extract_from_reaction
import pandas as pd
from collections import Counter
from tqdm import tqdm
import os
import argparse

def extract_from_rxn_smiles(rxn_smiles):
    reac, prod = rxn_smiles.split(">>")
    reaction_dict = {"reactants":reac, "products": prod, "_id":0, }
    template_dict = extract_from_reaction(reaction_dict)
    return template_dict.get("reaction_smarts")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--reaction_type',type=str)
    args = parser.parse_args()
    reaction_type = args.reaction_type

    print(reaction_type)
    rdchiral_templates = []
    data_dir = "../data/{}/".format(reaction_type)
    data_path = "../data/{}/MappingResult_{}.txt".format(reaction_type, reaction_type)
    with open(data_path, "r") as f:
        data = f.readlines()
    for line in tqdm(data):
        rxn_smiles = line.strip().split("\t")[0]
        template = extract_from_rxn_smiles(rxn_smiles)
        if not template:
            continue
        rdchiral_templates.append(template)
    template_count = dict(Counter(rdchiral_templates))
    for key, values in template_count.items():
        template_count[key] = [values]
    template_df  = pd.DataFrame.from_dict(template_count).transpose().reset_index().rename(columns={0:"count", "index":"template"})
    template_df = template_df.sort_values(by="count",ascending=False)
    template_df.to_csv(os.path.join(data_dir, "rdchiral_templates.csv"),index=False)
