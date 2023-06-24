"""
Use rxnmapper to give atom-mapping for the original Reaxys fiels.
"""

import os
import pandas as pd
import argparse
from rdkit import Chem
from rxnmapper import RXNMapper



def ValidRxnSmiles(rxn_smiles, reduce_radical = False):
    try:
        reac, prod = rxn_smiles.split('>>')
    except:
        return False
    if (reac == '') or (prod == ''):
        return False
    reac_mol = Chem.MolFromSmiles(reac)
    if reac_mol == None:
        return False
    for atom in reac_mol.GetAtoms():
        atom.SetIsotope(0)
        if reduce_radical:
            if atom.GetNumRadicalElectrons() != 0:
                atom.SetNumRadicalElectrons(0)
    prod_mol = Chem.MolFromSmiles(prod)
    if prod_mol == None:
        return False
    for atom in prod_mol.GetAtoms():
        atom.SetIsotope(0)
        if reduce_radical:
            if atom.GetNumRadicalElectrons() != 0:
                atom.SetNumRadicalElectrons(0)
    rxn_smiles = Chem.MolToSmiles(reac_mol) + '>>' + Chem.MolToSmiles(prod_mol)
    return rxn_smiles



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir',type=str,
                        default='../rxn_yield_context/Data_From_Reaxys_Original/FischerIndoleSynthesis')
    parser.add_argument('--output_file',type=str,
                        default='./data/FischerIndoleSynthesis/MappingResult_FischerIndoleSynthesis.txt')
    parser.add_argument('--reduce_radical', default=False, action='store_true')
    args = parser.parse_args()
    
    input_file_path = args.input_dir
    output_file_path = args.output_file
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
    
    
    input_files  = [os.path.join(input_file_path, x) for x in os.listdir(input_file_path) if x.endswith('.xlsx')]
    input_file = input_files[0]
    
    data = pd.read_excel(input_file, engine='openpyxl')
    keep_index = ['Reaction ID', 'Reaction']
    drop_index = list(data.columns)
    for index in list(data.columns):
        if index in keep_index:
            drop_index.remove(index)
    
    data = data.drop(drop_index,axis=1)
    data = data[:-3] # copyright or something 
    
    '''
    Concatenate all the data in the input_file_path directory. 
    '''
    for input_file in input_files[1:]:
        df1 = pd.read_excel(input_file, engine='openpyxl')
        df1 = df1.drop(drop_index, axis=1)
        df1 = df1[:-3]
        data = data.append(df1, ignore_index=True)
    
    print('Raw data information: ')
    data.info()
    print('-'*50 + '\n')
    
    '''
    Outputs atom-mapped rxn smiles.
    '''
    rxn_mapper = RXNMapper()
    reaxys_id_record = []
    processed_data = []
    for i in range(len(data)):
        reaxys_id = data.iloc[i]['Reaction ID']
        if reaxys_id in reaxys_id_record:
            continue
        else:
            reaxys_id_record.append(reaxys_id)
        rxn_smiles = data.iloc[i]['Reaction']
        new_rxn_smiles = ValidRxnSmiles(rxn_smiles, args.reduce_radical)
        if not new_rxn_smiles: continue
        
        try:
            results = rxn_mapper.get_attention_guided_atom_maps([new_rxn_smiles])
            mapped_rxn_smiles = results[0]['mapped_rxn']
            processed_data.append(mapped_rxn_smiles + '\t' + reaxys_id + '\n')
        except:
            continue
    
    f = open(output_file_path, 'w')
    f.writelines(processed_data)
    f.close()
        