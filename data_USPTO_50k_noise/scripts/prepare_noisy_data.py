""" This script aims to prepare the noisy dataset by (1) removing one reactant, (2) randomly changing the product atom and (3) randomly changing the atom mapping 
in reaction SMILES of USPTO-50k dataset.
`USPTO-50k golden dataset sourced from: https://raw.githubusercontent.com/connorcoley/retrosim/master/retrosim/data/data_processed.csv `
"""
import sys
sys.path.append("../../")
from rdkit import Chem
from autotemplate.run_utils import rdchiralRunText_modified
from tqdm import tqdm

import CGRtools
import random
import multiprocessing
import argparse
import os

def cgr_representation(rxn_smiles):
    r = CGRtools.smiles(rxn_smiles)
    r.clean2d()
    cgr = ~r
    return cgr.__str__()

def filter_invalid_outcomes(outcomes):
    valid_outcomes = []
    for smiles in outcomes:
        if Chem.MolFromSmiles(smiles):
            valid_outcomes.append(smiles)
    return valid_outcomes

def remove_one_reactant(rxn_smiles):
    reactants, product = rxn_smiles.split(">>")
    if "." not in reactants:
        return rxn_smiles
    else:
        reactants = reactants.split(".")
        reactant = random.choice(reactants)
        return reactant + ">>" + product

def random_change_product_atom(rxn_smiles):
    """May produce none list."""
    template_pools = ["[#6:1]>>[#6:1]-[#6]", # add a carbon atom 
                      "[#6:1]>>[#6:1]-[#8]", # add an oxygen atom 
                      "[#6:1]-[#6;D2:2]-[#6:3]>>[#6:1]-[#8;+0:4]-[#6:3]", # replace a carbon atom with an oxygen atom
                      "[#6:1]-[#6;D2:2]-[#6:3]>>[#6:1]-[#7;+0:4]-[#6:3]", # replace a carbon atom with an nitrogen atom
                      ]
    reactant, product = rxn_smiles.split(">>")
    t_index = list(range(len(template_pools)))
    random.shuffle(t_index)
    for id_ in t_index:
        random_products = rdchiralRunText_modified(template_pools[id_], product)
        random_products = filter_invalid_outcomes(random_products)
        if random_products:
            return reactant + ">>" + random.choice(random_products)
    else:
        print("This rxn smiles cannot be mutated:\n", rxn_smiles)
        return rxn_smiles

def random_change_atom_mapping(rxn_smiles):
    """Only needs to change the atom-mapping in the product. """
    reactants, product = rxn_smiles.split(">>")
    product_mol = Chem.MolFromSmiles(product)
    num_atom_dict = {}
    for atom in product_mol.GetAtoms():
        num_atom_dict.update({atom.GetAtomMapNum(): atom})

    if product_mol.GetNumAtoms() < 2:
        return rxn_smiles
    cgr_1 = cgr_representation(rxn_smiles)
    cgr_2 = cgr_1
    i = 0
    while cgr_1 == cgr_2:
        map_1, map_2 = random.sample(list(num_atom_dict.keys()), 2)
        atom_1 = num_atom_dict[map_1]
        atom_2 = num_atom_dict[map_2]
        if atom_1.GetAtomicNum() != atom_2.GetAtomicNum(): # Remember that the two changed atoms must be the same.
            continue
        num_atom_dict[map_1].SetAtomMapNum(map_2)
        num_atom_dict[map_2].SetAtomMapNum(map_1)
        mutated_product_smiles = Chem.MolToSmiles(product_mol)
        new_rxn_smiles = reactants + ">>" + mutated_product_smiles
        cgr_2 = cgr_representation(new_rxn_smiles)
        i += 1
        if i == 100:
            return rxn_smiles

    return new_rxn_smiles

def mp_wrapper(lines, task, ncpus):
    lines = [line.split("\t") for line in lines]
    rxn_smiles_list, reaction_id_list = list(zip(*lines))
    if task == "remove_one_reactant":
        subfunction = remove_one_reactant
    elif task == "random_change_product_atom":
        subfunction = random_change_product_atom
    elif task == "random_change_atom_mapping":
        subfunction = random_change_atom_mapping
    
    with multiprocessing.Pool(processes=ncpus) as pool:
        modified_rxn_smiles_list = pool.map(subfunction, rxn_smiles_list)

    return [modified_rxn_smiles + "\t" + reaction_id for modified_rxn_smiles, reaction_id in zip(modified_rxn_smiles_list, reaction_id_list)]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ncpus', type=int, default=16)
    args = parser.parse_args()
    random.seed(42)

    golden_data_path = "../../data_USPTO_50k_coley/USPTO_50k_coley.txt"
    output_data_path = "../"
    fractions = [3, 5, 10, 15, 20]

    for fraction in tqdm(fractions):
        with open(golden_data_path, "r") as f:
            input_lines = f.readlines()
        num_data = len(input_lines)
        num_noisy = num_data * fraction // 100
        index_list = list(range(num_data))

        # First task, remove one reactant:
        print("First task")
        effective_index_list = [i for i in index_list if "." in input_lines[i].split(">>")[0]]
        remove_index = sorted(random.sample(effective_index_list, num_noisy))
        remove_lines = [input_lines[i] for i in remove_index]
        results = mp_wrapper(remove_lines, "remove_one_reactant", args.ncpus)
        for remove_id, result in zip(remove_index, results):
            input_lines[remove_id] = result

        effective_index_list = [x for x in index_list if x not in remove_index] # remained index after the remove index

        # Second task, randomly change a product atom:
        print("Second task")
        change_atom_index = sorted(random.sample(effective_index_list, num_noisy))
        change_atom_lines = [input_lines[i] for i in change_atom_index]
        results = mp_wrapper(change_atom_lines, "random_change_product_atom", args.ncpus)
        for change_atom_id, result in zip(change_atom_index, results):
            input_lines[change_atom_id] = result
        
        effective_index_list = [x for x in index_list if x not in remove_index + change_atom_index] # remained index after the change atom index

        # Third task, randomly change the atom-mapping in product:
        print("Third task")
        change_mapping_index = sorted(random.sample(effective_index_list, num_noisy))
        change_mapping_lines = [input_lines[i] for i in change_mapping_index]
        results = mp_wrapper(change_atom_lines, "random_change_atom_mapping", args.ncpus)
        for change_mapping_id, result in zip(change_mapping_index, results):
            input_lines[change_mapping_id] = result

        # Ouput the file:
        os.makedirs(f"../USPTO_50k_{fraction}perc_noise", exist_ok=True)
        with open(f"../USPTO_50k_{fraction}perc_noise/USPTO_50k_{fraction}perc_noise.txt", "w") as g:
            g.writelines(input_lines)