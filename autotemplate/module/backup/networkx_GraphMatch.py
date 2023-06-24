# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 16:20:38 2021

@author: Lung-Yi

https://stackoverflow.com/questions/54476702/isomorphism-in-networkx-with-edge-attributes
https://networkx.org/documentation/stable/reference/algorithms/isomorphism.vf2.html
"""
import networkx as nx
import re
from networkx.algorithms import isomorphism
from rdkit import Chem

def canon_remap_local(smiles):
    # copy mol before changing it
    mol = Chem.MolFromSmiles(smiles)
    if mol == None: return None
    # remove all atom numbers beforehand, as the CanonicalRankAtoms function
    # takes into account existing mapping numbers while doing it's mapping (!)
    for atom in mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            atom.ClearProp('molAtomMapNumber')
    return Chem.MolToSmiles(mol,isomericSmiles=False)

def get_atoms_attributes(rdkit_molecule):
    attributes = []
    for a in rdkit_molecule.GetAtoms():
        attributes.append((a.GetAtomMapNum(), {'Symbol':a.GetSymbol(), 'Charge':a.GetFormalCharge()}))
    return attributes

def get_bonds_attributes(rdkit_molecule):
    attributes = []
    for b in rdkit_molecule.GetBonds():
        start = b.GetBeginAtom().GetAtomMapNum()
        end = b.GetEndAtom().GetAtomMapNum()
        b_type = b.GetBondTypeAsDouble()
        attributes.append((start, end, b_type))
    return attributes

def topology_from_rdkit(rdkit_molecule):
    topology = nx.Graph()
    topology.add_nodes_from(get_atoms_attributes(rdkit_molecule))
    topology.add_weighted_edges_from(get_bonds_attributes(rdkit_molecule))
    return topology


def mapping_for_gold_multiple_smiles(gold_smiles, outcome_smiles):
    # # !TODO: change this function to process two pairs of molecules.
    # gold_smiles = '[cH:3]1[c:13]([OH:25])[cH:14][cH:15][c:16]([CH:17]=[O:8])[c:29]1[OH:30].NNCc1ccc(O)cc1O'
    # # gold_smiles = 'NNCc1ccc(O)cc1O'
    # outcome_smiles ='[OH:12][c:13]1[cH:14][cH:15][c:16]([CH:17]=[O:21])[c:29]([OH:30])[cH:31]1.[NH2:18][NH:19][CH2:20][c:21]1[cH:22][cH:23][c:24]([OH:25])[cH:26][c:27]1[OH:28]'
    # # outcome_smiles = '[NH2:18][NH:19][CH2:20][c:21]1[cH:22][cH:23][c:24]([OH:25])[cH:26][c:27]1[OH:28]'
    # gold_smiles = 'O=[C:2]([CH2:3][c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1)[NH2:10].[NH2:12][C:13](=[O:14])[CH2:15][c:16]1[cH:17][cH:18][cH:19][cH:20][cH:21]1.[O:1]=[CH:11][c:22]1[cH:23][cH:24][c:25]([Br:26])[cH:27][cH:28]1'
    # outcome_smiles = '[CH:11]([c:22]1[cH:23][cH:24][c:25]([Br:26])[cH:27][cH:28]1)=[O:29].[NH2:12][C:13](=[O:14])[CH2:15][c:16]1[cH:17][cH:18][cH:19][cH:20][cH:21]1.[O:1]=[C:2]([CH2:3][c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1)[NH2:10]'
    if ('.' in gold_smiles) and ('.' in outcome_smiles):
        gold_split = gold_smiles.split('.')
        outcome_split = outcome_smiles.split('.')
        mapping_dict = dict(zip([smiles for smiles in gold_split], [{'canon': canon_remap_local(smiles)} for smiles in gold_split]))
        for smiles in outcome_split:
            outcome_canon = canon_remap_local(smiles)
            for key in mapping_dict:
                if (mapping_dict[key]['canon'] == outcome_canon) and (mapping_dict[key].get('outcome') == None):
                    mapping_dict[key].update({'outcome': smiles})
                    break
        # [mapping_dict[canon_remap(smiles)].update({'outcome': smiles}) for smiles in outcome_split]
        final_smiles_list = []
        for key, value in mapping_dict.items():
            smiles = mapping_for_gold_smiles(key, value.get('outcome'), show_info = False)
            if smiles:
                final_smiles_list.append(smiles)
        return '.'.join(final_smiles_list)
    else:
        return mapping_for_gold_smiles(gold_smiles, outcome_smiles)

def mapping_for_gold_smiles(gold_smiles, outcome_smiles, show_info = False):
    """
    Parameters
    ----------
    gold_smiles : str
        gold_smiles is Reaxys smiles with stereo and chiral information.
    outcome_smiles : str
        outcome_smiles is the reactant or product after template-mapped,
        and it does not contain stereo and chiral information.

    Returns
    -------
    gold_smiles : str

    """
    
    if (gold_smiles == None) or (outcome_smiles == None):
        return
    gold_mol = Chem.MolFromSmiles(gold_smiles)
    outcome_mol = Chem.MolFromSmiles(outcome_smiles)
    # Chem.SanitizeMol(gold_mol)
    # Chem.SanitizeMol(outcome_mol)
    for a in gold_mol.GetAtoms():
        a.SetAtomMapNum(a.GetIdx())
    gold_Graph = topology_from_rdkit(gold_mol)
    outcome_Graph = topology_from_rdkit(outcome_mol)
    GM = isomorphism.GraphMatcher(gold_Graph, outcome_Graph,
                              node_match=lambda n1,n2:n1['Symbol']==n2['Symbol'],
                              edge_match=lambda e1,e2:e1['weight']==e2['weight'])
    if (not GM.is_isomorphic()) : # GM.is_isomorphic() function can only consider one pair of molecules at the same time.
        if show_info:
            print('-'*60)
            print('The reaction outcome and gold answer do not match graph isomorphism!')
            print('gold_smiles:')
            print(Chem.MolToSmiles(gold_mol))
            print('outcome_smiles:')
            print(outcome_smiles)
        return
    else:
        mapping = GM.mapping
        for a in gold_mol.GetAtoms():
            a.SetAtomMapNum(mapping[a.GetAtomMapNum()])
        # print(Chem.MolToSmiles(gold_mol))
        return Chem.MolToSmiles(gold_mol)

def compare_two_graphs(first, second):
    GM_reac = isomorphism.GraphMatcher(first['reactant'], second['reactant'], 
                                       node_match=isomorphism.categorical_node_match(['Symbol', 'Charge'], ['C',0]),
                                       edge_match=lambda e1,e2:e1['weight']==e2['weight'])
    GM_prod = isomorphism.GraphMatcher(first['product'], second['product'], 
                                       node_match=isomorphism.categorical_node_match(['Symbol', 'Charge'], ['C',0]),
                                       edge_match=lambda e1,e2:e1['weight']==e2['weight'])
    return GM_reac.is_isomorphic() and GM_prod.is_isomorphic()

def SmartsMolWithFormalCharge(smarts):
    """
    [O+]-[#6:3]=[#7;+1:4]=[#7;-2:5]
    """
    mol = Chem.MolFromSmarts(smarts)
    # smarts_convert = ']'.join([x for x in smarts.split('[') if (':' in x) and (x != ':')])
    # smarts_convert = [x for x in smarts_convert.split(']') if (':' in x) and (x != ':')]
    smarts_convert = re.findall('\[([^]]*)\]', smarts)
    smarts_convert = [x for x in smarts_convert if ':' in x]
    formalcharge_dict = dict()
    for symbol in smarts_convert:
        charge = re.search('([-+]+[1-9]?)', symbol)
        if charge:
            mapnum = int(re.search('\:([0-9]+)', symbol).group()[1:])
            formalcharge_dict.update({mapnum: int(charge.group())})
    
    for atom in mol.GetAtoms():
        mapnum = int(atom.GetAtomMapNum())
        if mapnum and (mapnum in formalcharge_dict.keys()):
            atom.SetFormalCharge(formalcharge_dict[mapnum])
    
    return mol

def convert_template2graphs(template):
    reac, prod = template.split('>>')
    reac_graph = topology_from_rdkit(SmartsMolWithFormalCharge(reac))
    prod_graph = topology_from_rdkit(SmartsMolWithFormalCharge(prod))
    return {'reactant': reac_graph, 'product': prod_graph}

def find_unique_templates(all_templates):
    unique_templates = [all_templates[0]]
    Graphs_Reaction_Change = [convert_template2graphs(all_templates[0])]
    for template in all_templates:
        # template = all_templates[1]
        new_graphs = convert_template2graphs(template)
        for old_graphs in Graphs_Reaction_Change:
            if compare_two_graphs(old_graphs, new_graphs): # duplicate reaction template in records
                break
        else:
            unique_templates.append(template)
            Graphs_Reaction_Change.append(new_graphs)
    return unique_templates

def find_unique_templates_dict(templates_dict):
    unique_templates = dict()
    changed_records = dict()
    # Graphs_Reaction_Change = [convert_template2graphs(all_templates[0])]
    for template, count_value in templates_dict.items():
        if not unique_templates:
            print('Init')
            information = convert_template2graphs(template)
            information.update({'count': count_value})
            unique_templates.update({template: information})
            
            changed_records.update({template: template})
        else:
            new_graphs = convert_template2graphs(template)
            for key, old_graphs in unique_templates.items():
                if compare_two_graphs(old_graphs, new_graphs): # duplicate reaction template in records
                    unique_templates[key]['count'] += count_value
                    
                    changed_records.update({template: key})
                    break
            else:
                information = convert_template2graphs(template)
                information.update({'count': count_value})
                unique_templates.update({template: information})
                
                changed_records.update({template: template})
    
    for key, value in unique_templates.items():
        unique_templates[key] = value['count']
    return unique_templates, changed_records
    
    
# if __name__ == '__main__':
#     pass
#     gold_smiles = '[CH3:1][CH2:2][O:3][C:4](=[O:5])/[CH:6]=[CH:7]/[C:8](=[O:9])[O:10][CH2:11][CH3:12].[CH:13]1=[CH:14][CH:15]=[CH:16][CH2:17]1'
#     outcome_smiles = '[CH3:1][CH2:2][O:3][C:4](=[O:5])[CH:6]=[CH:12][C:13](=[O:14])[O:15][CH2:16][CH3:17].[CH:7]1=[CH:8][CH:9]=[CH:10][CH2:11]1'
    
#     print('gold_smiles:')
#     print(gold_smiles)
#     print('outcome_smiles')
#     print(outcome_smiles)
#     print('gold_smiles after remapped:')
#     print(mapping_for_gold_smiles('[CH3:1][CH2:2][O:3][C:4](=[O:5])/[CH:6]=[CH:7]/[C:8](=[O:9])[O:10][CH2:11][CH3:12].[CH:13]1=[CH:14][CH:15]=[CH:16][CH2:17]1',
#                           '[CH3:1][CH2:2][O:3][C:4](=[O:5])[CH:6]=[CH:12][C:13](=[O:14])[O:15][CH2:16][CH3:17].[CH:7]1=[CH:8][CH:9]=[CH:10][CH2:11]1'))
