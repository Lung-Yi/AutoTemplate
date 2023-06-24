# -*- coding: utf-8 -*-
"""
https://stackoverflow.com/questions/54476702/isomorphism-in-networkx-with-edge-attributes
https://networkx.org/documentation/stable/reference/algorithms/isomorphism.vf2.html
"""
import networkx as nx
import re
from networkx.algorithms import isomorphism
from rdkit import Chem
import numpy as np

def clear_atom_map_rxn_smiles(rxn_smiles):
    reac, prod = rxn_smiles.split('>>')
    reac = canon_remap_local(reac, True)
    prod = canon_remap_local(prod, True)
    return reac + '>>' + prod

def canon_remap_local(smiles, iso = False):
    # copy mol before changing it
    mol = Chem.MolFromSmiles(smiles)
    if mol == None: return None
    # remove all atom numbers beforehand, as the CanonicalRankAtoms function
    # takes into account existing mapping numbers while doing it's mapping (!)
    for atom in mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            atom.ClearProp('molAtomMapNumber')
    return Chem.MolToSmiles(mol,isomericSmiles=iso)

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


""" For template canonicalizaiton """

def get_atoms_attributes_idx(rdkit_molecule):
    attributes = []
    for a in rdkit_molecule.GetAtoms():
        attributes.append((a.GetIdx(), {'Symbol':a.GetSymbol(), 'Charge':a.GetFormalCharge()}))
    return attributes

def get_bonds_attributes_idx(rdkit_molecule):
    attributes = []
    for b in rdkit_molecule.GetBonds():
        start = b.GetBeginAtom().GetIdx()
        end = b.GetEndAtom().GetIdx()
        b_type = b.GetBondTypeAsDouble()
        attributes.append((start, end, b_type))
    return attributes

def topology_from_rdkit_idx(rdkit_molecule):
    topology = nx.Graph()
    topology.add_nodes_from(get_atoms_attributes_idx(rdkit_molecule))
    topology.add_weighted_edges_from(get_bonds_attributes_idx(rdkit_molecule))
    return topology


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
    """

    Parameters
    ----------
    template : str
        A retro reaction template (reaction SMARTS).

    Returns
    -------
    dict
        Three variables: 'reactant'-> reactant graph; 'product'-> product graph; 'BCN'-> Bond change number of this reaction.

    """
    reac, prod = template.split('>>')
    reac_graph = topology_from_rdkit_idx(SmartsMolWithFormalCharge(reac))
    prod_graph = topology_from_rdkit_idx(SmartsMolWithFormalCharge(prod))
    BCN = countBCN(template)
    return {'reactant': reac_graph, 'product': prod_graph, 'BCN': BCN}

def find_unique_templates(all_templates):
    '''' Currently this funciton is not used. '''
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
    """ This function aims to:
        (1) merge the same reaction template by checking the Graph isomorphism.
        (2) remove the template that the product site has more atoms. (It is allowed that reactant site has more atoms.)
        """
    unique_templates = dict()
    changed_records = dict()
    sorted_template_values = sorted(list(templates_dict.items()), key=lambda x:x[1], reverse=True)
    # Graphs_Reaction_Change = [convert_template2graphs(all_templates[0])]
    for template, count_value in sorted_template_values:
        p, r = template.split('>>')
        p_atoms = Chem.MolFromSmarts(p).GetNumAtoms()
        r_atoms = Chem.MolFromSmarts(r).GetNumAtoms()
        if p_atoms > r_atoms:
            continue
        
        if not unique_templates:
            print('Init')
            information = convert_template2graphs(template)
            information.update({'count': count_value})
            unique_templates.update({template: information})
            
            changed_records.update({template: template})
        else:
            new_graphs = convert_template2graphs(template)
            if new_graphs['BCN'] == np.inf:
                continue
                
            for key, old_graphs in unique_templates.items():
                if compare_two_graphs(old_graphs, new_graphs): # duplicate reaction template in records
                    if old_graphs['BCN'] <= new_graphs['BCN']: # merge new template into old template
                        unique_templates[key]['count'] += count_value
                        changed_records.update({template: key})
                        break
                    elif old_graphs['count'] < count_value: # replace the old template with the new template because the template with less BCN is correct
                        new_information = convert_template2graphs(template)
                        new_information.update({'count': count_value + old_graphs['count']})
                        unique_templates.update({template: new_information})
                        del unique_templates[key]
                        changed_records.update({key: template})
                        break
            else:
                information = convert_template2graphs(template)
                information.update({'count': count_value})
                unique_templates.update({template: information})
                changed_records.update({template: template})
    
    for key, value in unique_templates.items():
        unique_templates[key] = value['count']
    return unique_templates, changed_records



def swap(i,j):
    return min(i,j), max(i,j)

def CheckTwoAtoms(atom_r, atom_p):
    atom_r_info = set()
    atom_p_info = set()
    for bond in atom_r.GetBonds():
        start = bond.GetBeginAtom().GetAtomMapNum()
        end = bond.GetEndAtom().GetAtomMapNum()
        start, end = swap(start,end)
        info = str(start) +'-' + str(end) +'-'+ str(bond.GetBondTypeAsDouble())
        atom_r_info.add(info)
    for bond in atom_p.GetBonds():
        start = bond.GetBeginAtom().GetAtomMapNum()
        end = bond.GetEndAtom().GetAtomMapNum()
        start, end = swap(start,end)
        info = str(start) +'-' + str(end) +'-'+ str(bond.GetBondTypeAsDouble())
        atom_p_info.add(info)
    if atom_r_info == atom_p_info:
        return
    else:
        answer = []
        p_check = list(atom_p_info)
        p_check = [i.split('-')[0:2] for i in p_check]
        for element in (atom_r_info ^ atom_p_info):
            if element in atom_p_info:
                answer.append(element)
            else:
                if element.split('-')[0:2] in p_check:
                    continue
                else:
                    e = element.split('-')
                    e[-1] = '0.0'
                    answer.append('-'.join(e))
        return answer

def GetAnswer(rxn_smiles):
    """ If reactant site has more unmapped atom, fill up it with unoccupied atom map number. """
    try:
        AllBondChange = []
        rea = rxn_smiles.split('>')[0]
        prod = rxn_smiles.split('>')[-1]
        r = Chem.MolFromSmarts(rea)
        p = Chem.MolFromSmarts(prod)
        map_num_pool = set([i+1 for i in range(r.GetNumAtoms() + p.GetNumAtoms())])
        for atom in r.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num in map_num_pool:
                map_num_pool.remove(map_num)
        for atom in p.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num in map_num_pool:
                map_num_pool.remove(map_num)
        
        r_info = {}
        p_info = {}
        for atom in r.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num == 0:
                map_num = map_num_pool.pop()
                atom.SetAtomMapNum(map_num)
            r_info.update({map_num: atom})
        for atom in p.GetAtoms():
            p_info.update({atom.GetAtomMapNum():atom})
        
        for key,value in p_info.items():
            atom_r = r_info[key]
            atom_p = value
            result = CheckTwoAtoms(atom_r, atom_p)
            if result:
                AllBondChange += result
        AllBondChange = list(set(AllBondChange))
        return ';'.join(AllBondChange)
    except:
        # print(rxn_smiles)
        return None

def CountNumBondChange(answer):
    if ';' not in answer:
        return 1
    else:
        return answer.count(';')+1

def countBCN(template):
    A, B = template.split('>>')
    forward_template = B + '>>' +A
    changed_list = GetAnswer(forward_template)
    if changed_list:
        return CountNumBondChange(changed_list)
    else:
        return np.inf
    
    
if __name__ == '__main__':
    pass
    gold_smiles = '[CH3:1][CH2:2][O:3][C:4](=[O:5])/[CH:6]=[CH:7]/[C:8](=[O:9])[O:10][CH2:11][CH3:12].[CH:13]1=[CH:14][CH:15]=[CH:16][CH2:17]1'
    outcome_smiles = '[CH3:1][CH2:2][O:3][C:4](=[O:5])[CH:6]=[CH:12][C:13](=[O:14])[O:15][CH2:16][CH3:17].[CH:7]1=[CH:8][CH:9]=[CH:10][CH2:11]1'
    
    print('gold_smiles:')
    print(gold_smiles)
    print('outcome_smiles')
    print(outcome_smiles)
    print('gold_smiles after remapped:')
    print(mapping_for_gold_smiles('[CH3:1][CH2:2][O:3][C:4](=[O:5])/[CH:6]=[CH:7]/[C:8](=[O:9])[O:10][CH2:11][CH3:12].[CH:13]1=[CH:14][CH:15]=[CH:16][CH2:17]1',
                          '[CH3:1][CH2:2][O:3][C:4](=[O:5])[CH:6]=[CH:12][C:13](=[O:14])[O:15][CH2:16][CH3:17].[CH:7]1=[CH:8][CH:9]=[CH:10][CH2:11]1'))
