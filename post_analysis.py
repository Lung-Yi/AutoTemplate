# -*- coding: utf-8 -*-

import pandas as pd
import os
from rdkit import Chem
from tqdm import tqdm
from autotemplate.run_utils import clearIsotope, RemoveReagent
from autotemplate.extract_utils import canon_remap
import matplotlib.pyplot as plt
import CGRtools
plt.rcParams["figure.dpi"] = 400
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def make_dict(data):
    index = dict()
    for line in data:
        rxn_smiles, rxn_id = line.strip().split('\t')
        index.update({rxn_id:rxn_smiles})
    return index

def CGR_smarts(reaction_smiles, display_figure=False):
    reaction_smiles_new = clearIsotope(reaction_smiles)
    r = CGRtools.smiles(reaction_smiles_new)
    r.clean2d()
    cgr = ~r
    if display_figure:
        display(cgr)
    return cgr.__str__()

def is_atommap_corrected(ori_rxn_smiles, proc_rxn_smiles):
    cgr_ori = CGR_smarts(ori_rxn_smiles)
    cgr_proc = CGR_smarts(proc_rxn_smiles)
    return not cgr_ori == cgr_proc


def is_curated(ori_rxn_smiles, proc_rxn_smiles):
    """ If the processed reactant SMILES contains the molecule that original reactant does not posses, 
    this reaction has been curated. Use unsorted list to check the elements. """
    ori_reac_smiles = ori_rxn_smiles.split('>>')[0]
    proc_reac_smiles = proc_rxn_smiles.split('>>')[0]
    ori_reac_list = [canon_remap(s) for s in ori_reac_smiles.split()]
    proc_reac_list = [canon_remap(s) for s in proc_reac_smiles.split()]
    for smi in proc_reac_list:
        if smi in ori_reac_list:
            ori_reac_list.remove(smi)
        else:
            return True
    else:
        return False

all_rxn_class = ["AdamsDecarboxylation",
                "Baylis-HillmanReaction",
                "Buchwald-HartwigCross-Coupling",
                "Chan_LamCoupling",
                "DielsAlder",
                "FischerIndoleSynthesis",
                "Friedel-CraftsAcylation",
                "Friedel-CraftsAlkylation",
                "GrignardReaction",
                "HiyamaCoupling",
                "HuisgenCycloaddition",
                "Hydrogenation",
                "Kabachnik-FieldsReaction",
                "KumadaCoupling",
                "MannichReaction",
                "NegishiCoupling",
                "PausonKhandReaction",
                "ReductiveAmination",
                "SuzukiCoupling",
                "WittigReaction",
                ]

records = []

for k, rxn_class in enumerate(all_rxn_class):
    print('Current rxn: {}, number {}'.format(rxn_class, k+1))
    data_dir = './data/{}'.format(rxn_class)
    unprocessed_path = os.path.join(data_dir, 'MappingResult_{}.txt'.format(rxn_class))
    processed_path = os.path.join(data_dir, 'MappingResult_{}.txt.processed'.format(rxn_class))
    failed_path = os.path.join(data_dir, 'MappingResult_{}.txt.failed'.format(rxn_class))
    
    with open(unprocessed_path, 'r') as f:
        unprocessed = f.readlines()
    with open(processed_path, 'r') as f:
        processed = f.readlines()
    with open(failed_path, 'r') as f:
        failed = f.readlines()
        
    total = len(unprocessed)
    same = 0
    curated = 0
    mapping_curated = 0
    failed = len(failed)
    unprocessed = make_dict(unprocessed)
    processed = make_dict(processed)
    
    for rxn_id, rxn_smiles_1 in tqdm(processed.items()):
        rxn_smiles_2 = unprocessed[rxn_id]
        rxn_smiles_2 = RemoveReagent(rxn_smiles_2)
        try:
            smarts_1 = CGR_smarts(rxn_smiles_1)
            smarts_2 = CGR_smarts(rxn_smiles_2)
        except:
            print(rxn_id)
            print(rxn_smiles_1)
            print(rxn_smiles_2)
            mapping_curated += 1
            continue
        if is_curated(rxn_smiles_1, rxn_smiles_2):
            curated += 1
        elif is_atommap_corrected(rxn_smiles_1, rxn_smiles_2):
            mapping_curated += 1
        else:
            same += 1
    
    records.append([same/total, mapping_curated/total, curated/total, failed/total])

total_data = pd.DataFrame(records, columns = ['no change', 'mapping curated' ,'reactant curated', 'removed'], index = all_rxn_class)
total_data = total_data.iloc[::-1]
total_data.to_csv('analyze_results.csv')
plot = total_data.plot(kind='barh',stacked=True)

fig = plot.get_figure()
fig.savefig("output.svg", format="svg")

