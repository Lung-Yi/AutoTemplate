# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 23:33:15 2020

@author: Lung-Yi
"""
from __future__ import print_function
import sys 
import os
import re
import copy
import math
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir

from rdchiral.utils import vprint, PLEVEL, atoms_are_different
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.chiral import template_atom_could_have_been_tetra, atom_chirality_matches#,copy_chirality,\
    
from rdchiral.clean import canonicalize_outcome_smiles, combine_enantiomers_into_racemic
from rdchiral.bonds import BondDirOpposite, restore_bond_stereo_to_sp2_atom

def canon_remap(smiles, return_NumAtom=False):
    # copy mol before changing it
    mol = Chem.MolFromSmiles(smiles)
    if mol == None: return None
    # remove all atom numbers beforehand, as the CanonicalRankAtoms function
    # takes into account existing mapping numbers while doing it's mapping (!)
    for atom in mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            atom.ClearProp('molAtomMapNumber')
    if return_NumAtom:
        return Chem.MolToSmiles(mol,isomericSmiles=False), mol.GetNumAtoms()
    return Chem.MolToSmiles(mol,isomericSmiles=False)

def CalculateNumHs(atom):
    # atom_valence_dict = {"C":4, "N":3, "O":2, "S":(2,4,6), "P":(3,5),"Si":4, "Br":1, "Mg":2}
    valence = atom_valence_dict[atom.GetSymbol()]
    return int(valence - sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]) - abs(atom.GetFormalCharge()))

# def UpdateAtomInfo(atom):
#     # atom_valence_dict = {"C":4, "N":3, "O":2, "S":(2,4,6), "P":(3,5),"Si":4, "Br":1, "Mg":2}
#     if atom.GetTotalValence() != atom_valence_dict[atom.GetSymbol()]:
#         atom.SetNumExplicitHs(CalculateNumHs(atom))
#         atom.UpdatePropertyCache()
#     return

def GetAtomWithAtomMapNum(mol, mapnum):
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == mapnum:
            return atom
    return

def ReassignMapping(smiles, return_NumAtom=False, iso=False):
    m = Chem.MolFromSmiles(smiles)
    i = 1
    for a in m.GetAtoms():
        if a.GetIsotope():
            a.SetIsotope(0)
            if a.GetSymbol() != 'H':
                a.SetAtomMapNum(i)
                i += 1  
        else:
            a.SetAtomMapNum(i)
            i += 1
    if return_NumAtom:
        return Chem.MolToSmiles(m,isomericSmiles=iso), m.GetNumAtoms()
    return Chem.MolToSmiles(m,isomericSmiles=iso)

atom_valence_dict = {"C":4, "N":3, "O":2, "S":(2,4,6), "P":(3,5),"Si":4, "Br":1, "Cl":1, "I":1, "F":1, "Mg":2}

"""
start testing
append_reagent function only supports one atom now
"""
def rdchiralRunText_modified(reaction_smarts, reactant_smiles, append_reagent=False):
    reactant_mapnum = set([int(atom.GetProp('molAtomMapNumber')) for atom in Chem.MolFromSmiles(reactant_smiles).GetAtoms()])
    rxn = rdchiralReaction(reaction_smarts)
    rxn.reset()
    if append_reagent:
        need_reagent = True
        if rxn.template_r.GetNumAtoms() == rxn.template_p.GetNumAtoms():
            need_reagent = False
        else:
            reagents_list = []
    reactants = rdchiralReactants(reactant_smiles)
    # Run naive RDKit on ACHIRAL version of molecules
    outcomes = rxn.rxn.RunReactants((reactants.reactants_achiral,))
    # mol = Chem.MolFromSmiles(reactant_smiles)
    # Chem.rdmolops.RemoveStereochemistry(mol)
    # outcomes = rxn.rxn.RunReactants((mol,)) # for tesret
    smiles_list = []
    for outcome in outcomes:
        ###############################################################################
        # Look for new atoms in products that were not in 
        # reactants (e.g., LGs for a retro reaction)
        changed_index = []
        reagent_smiles = []
        unmapped = sum([z.GetNumAtoms() for z in outcome])
        for m in outcome:
            try:
                for a in m.GetAtoms():
                    # if not a.IsInRing(): 
                    #     a.SetIsAromatic(False)
                    #     for b in a.GetBonds():
                    #         b.SetIsAromatic(False)
                    # Assign map number to outcome based on react_atom_idx
                    if a.HasProp('react_atom_idx'):
                        num = reactants.idx_to_mapnum(int(a.GetProp('react_atom_idx')))
                        if a.GetAtomMapNum() != num: 
                            changed_index.append(num)
                            a.UpdatePropertyCache(strict=False)
                        a.SetAtomMapNum(num)
            except:
                continue
            try:
                Chem.Kekulize(m)
            except Exception as e:
                # print(e)
                continue
            for a in m.GetAtoms():
                if (int(a.GetAtomMapNum()) in changed_index) and (not a.HasProp('_QueryHCount')): # check type
                    try:
                        if type(atom_valence_dict[a.GetSymbol()]) == int:
                            if a.GetTotalValence() != atom_valence_dict[a.GetSymbol()]:
                                a.SetNumExplicitHs(CalculateNumHs(a))
                                a.UpdatePropertyCache(strict=True)
                        else:
                            for valence_number in atom_valence_dict[a.GetSymbol()]:
                                try:
                                    a.SetNumExplicitHs(int(valence_number - sum([bond.GetBondTypeAsDouble() for bond in a.GetBonds()]) - abs(a.GetFormalCharge()) ))
                                    a.UpdatePropertyCache(strict=True)
                                    break
                                except:
                                    continue
                    except Exception as e:
                        # print(e)
                        pass
            try:
                Chem.SanitizeMol(m)
            except:
                continue
            
        #######################################################################################
        # Convert product(s) to single product so that all 
        # reactions can be treated as pseudo-intramolecular
        # But! check for ring openings mistakenly split into multiple
        # This can be diagnosed by duplicate map numbers (i.e., SMILES)
        mapnums = [a.GetAtomMapNum() for m in outcome for a in m.GetAtoms() if a.GetAtomMapNum()]
        if len(mapnums) != len(set(mapnums)): # duplicate?
            try:
                if PLEVEL >= 1: print('Found duplicate mapnums in product - need to stitch')
                # need to do a fancy merge
                merged_mol = Chem.RWMol(outcome[0])
                merged_map_to_id = {a.GetAtomMapNum(): a.GetIdx() for a in outcome[0].GetAtoms() if a.GetAtomMapNum()}
                for j in range(1, len(outcome)):
                    new_mol = outcome[j]#
                    for a in new_mol.GetAtoms():
                        if a.GetAtomMapNum() not in merged_map_to_id:
                            merged_map_to_id[a.GetAtomMapNum()] = merged_mol.AddAtom(a)
                    for b in new_mol.GetBonds():
                        bi = b.GetBeginAtom().GetAtomMapNum()
                        bj = b.GetEndAtom().GetAtomMapNum()
                        if PLEVEL >= 10: print('stitching bond between {} and {} in stich has chirality {}, {}'.format(
                            bi, bj, b.GetStereo(), b.GetBondDir()
                        ))
                        if not merged_mol.GetBondBetweenAtoms(
                                merged_map_to_id[bi], merged_map_to_id[bj]):
                            merged_mol.AddBond(merged_map_to_id[bi],
                                merged_map_to_id[bj], b.GetBondType())
                            merged_mol.GetBondBetweenAtoms(
                                merged_map_to_id[bi], merged_map_to_id[bj]
                            ).SetStereo(b.GetStereo())
                            merged_mol.GetBondBetweenAtoms(
                                merged_map_to_id[bi], merged_map_to_id[bj]
                            ).SetBondDir(b.GetBondDir())
                outcome = merged_mol.GetMol()
            except:
                continue
        else:
            new_outcome = outcome[0]
            for j in range(1, len(outcome)):
                new_outcome = AllChem.CombineMols(new_outcome, outcome[j])
            outcome = new_outcome
        
        # Assign atom map num after stitching to those atoms that are reagent atom:
        # "reactant_mapnum" is the set containing the atom mapping number that have been used in reactants
        unmapped_mapnum = set([ i+1 for i in range(0, 299)])
        unmapped_mapnum = unmapped_mapnum.difference(reactant_mapnum)
        # outcome_mapnum = set([atom.GetProp('molAtomMapNumber') for atom in Chem.MolFromSmiles(outcome).GetAtoms()])
        # unmapped_mapnum = unmapped_mapnum.difference(outcome_mapnum)
        
        for a in outcome.GetAtoms():
            if (not a.HasProp('react_atom_idx')) and (not a.GetAtomMapNum()):
                unmapped = unmapped_mapnum.pop()
                a.SetAtomMapNum(unmapped)
                if append_reagent:
                    reagent_atom = '[{}:{}]'.format(a.GetSymbol(), unmapped)
                    reagent_mol = Chem.MolFromSmiles(reagent_atom)
                    reagent_mol.GetAtoms()[0].SetNumExplicitHs(CalculateNumHs(reagent_mol.GetAtoms()[0]))
                    reagent_smiles.append(Chem.MolToSmiles(reagent_mol))
                # unmapped -= 1
                
        #######################################################################################
        # Update molecule because maybe bonds change.
        try:
            outcome.UpdatePropertyCache(strict=True)
        except:
            pass
        #######################################################################################
        smiles = Chem.MolToSmiles(outcome)
        if append_reagent:
            if need_reagent:
                reagent_smiles = '.'.join(reagent_smiles)
                reagents_list.append(reagent_smiles)
        smiles_list.append(smiles.rstrip("."))
    if append_reagent:
        if not need_reagent:
            reagents_list = ['']*len(smiles_list)
        return smiles_list, reagents_list
    return list(set(smiles_list))

def RemoveReagent(rxn_smiles, select_major_product = False):
    r, p = rxn_smiles.split('>>')
    r = r.split('.')
    # if '.' not in p:
    #     p = Chem.MolFromSmiles(p)
    #     p_map_total = [x.GetAtomMapNum() for x in p.GetAtoms()]
    # else:
    p = p.split('.')
    can_r = [Chem.CanonSmiles(smi) for smi in r]
    can_p = [Chem.CanonSmiles(smi) for smi in p]
    inter = set(can_r) & set(can_p)
    if inter:
        for smi in inter:
            can_p.remove(smi)
            can_r.remove(smi)
    r = can_r
    p = can_p
    if (p == []) or (r == []):
        print("This reaction has no change:")
        print(rxn_smiles)
        return rxn_smiles
    r = [Chem.MolFromSmiles(smi) for smi in r]
    p = [Chem.MolFromSmiles(smi) for smi in p]
    r_map_total = [atom.GetAtomMapNum() for mol in r for atom in mol.GetAtoms()]
    p_map_total = [atom.GetAtomMapNum() for mol in p for atom in mol.GetAtoms()]
    r2 = r.copy()   
    p2 = p.copy()
    for m in r2:
        r_map = [a.GetAtomMapNum() for a in m.GetAtoms()]
        if set(r_map) & set(p_map_total) == set(): r.remove(m)
    
    for m in p2:
        p_map = [a.GetAtomMapNum() for a in m.GetAtoms()]
        if set(r_map_total) & set(p_map) == set(): p.remove(m)
    
    if select_major_product:
        p = sorted(p, key= lambda x: x.GetNumAtoms(),reverse=True)[0] # select the major product
        return '.'.join([Chem.MolToSmiles(m) for m in r]) +'>>'+ Chem.MolToSmiles(p)
    else:
        return '.'.join([Chem.MolToSmiles(m) for m in r]) +'>>'+ '.'.join([Chem.MolToSmiles(x) for x in p])
        
# def Get_chirality(atom):
#     """Input an atom, and return the map number and its chirality info.
#        None means this atom is not chiral center."""
#     # Not possible to be a tetrahedral center anymore?
#     if atom.GetDegree() < 3:
#         return 
#     if atom.GetDegree() == 3 and \
#             any(b.GetBondType() != BondType.SINGLE for b in atom.GetBonds()):
#         return
    
#     return atom.GetAtomMapNum(), atom.GetChiralTag()


def copy_chirality_modify(a_src, a_new):
    """append chiral info to new correspnding atom """
    # Not possible to be a tetrahedral center anymore?
    if a_new.GetDegree() < 3:
        return 
    if a_new.GetDegree() == 3 and \
            any(b.GetBondType() != BondType.SINGLE for b in a_new.GetBonds()):
        return

    if PLEVEL >= 3: print('For mapnum {}, copying src {} chirality tag to new'.format(
        a_src.GetAtomMapNum(), a_src.GetChiralTag()))
    a_new.SetChiralTag(a_src.GetChiralTag())

def copy_stereo(b_src, b_new):
    """Notice: https://github.com/rdkit/rdkit/issues/2404 
       Need to append SetStereoAtoms(end_atom_idx, begin_atom_idx)"""
    if b_new.GetBondTypeAsDouble() == Chem.rdchem.BondStereo.STEREONONE:
        return
    b_new.SetStereoAtoms(b_new.GetEndAtomIdx(), b_new.GetBeginAtomIdx())
    b_new.SetStereo(b_src.GetStereo())
    return

def move_info(gold_smiles, new_smiles):
    g_count = gold_smiles.count('.')
    n_count = new_smiles.count('.')
    assert g_count == n_count
    if g_count == 0:
        gold_mol = Chem.MolFromSmiles(gold_smiles)
        new_mol = Chem.MolFromSmiles(new_smiles)
        append_ChiralStereo_info_for_mol(gold_mol,new_mol)
    else:
        gold_smiles = gold_smiles.split('.')
        new_smiles = new_smiles.split('.')
        gold_smiles = sorted(gold_smiles,key= lambda s:canon_remap(s))
        new_smiles = sorted(new_smiles,key= lambda s:canon_remap(s))
        for i in range(g_count+1):
            m1 = Chem.MolFromSmiles(gold_smiles[i])
            m2 = Chem.MolFromSmiles(new_smiles[i])
            append_ChiralStereo_info_for_mol(m1, m2)
            if i == 0:
                new_mol = m2
            else:
                new_mol = AllChem.CombineMols(new_mol, m2)
    return Chem.MolToSmiles(new_mol)
            
    

def append_ChiralStereo_info_for_mol(gold_mol, new_mol):
    """
    Parameters
    ----------
    gold_mol : rdkit.MolObject
        gold_mol has the chiral and stereo information.
    new_mol : rdkit.MolObject
        new_mol is the reactant of product after template-mapping, so it
        does not have chiral and stereo information in itself.

    Returns None
    -------
    None.

    """
    assert gold_mol.GetNumAtoms() == new_mol.GetNumAtoms()
    Chem.AssignStereochemistry(new_mol, force=True, cleanIt=True)
    gold_atoms = [atom for atom in gold_mol.GetAtoms()]
    new_atoms = [atom for atom in new_mol.GetAtoms()]
    for i in range(gold_mol.GetNumAtoms()):
        gold_atom = gold_atoms[i]
        new_atom = new_atoms[i]
        assert not atoms_are_different_2nd(gold_atom, new_atom)
        copy_chirality_modify(gold_atom, new_atom)
        gold_bonds = [bond for bond in gold_atom.GetBonds()]
        new_bonds = [bond for bond in new_atom.GetBonds()]
        gold_bonds = sorted(gold_bonds, key= lambda b:bond_to_label_2nd(b))
        new_bonds =  sorted(new_bonds, key= lambda b:bond_to_label_2nd(b))
        for j in range(len(gold_bonds)):
            copy_stereo(gold_bonds[j], new_bonds[j])
    #################Check the results#################
    # gold_smiles = Chem.MolToSmiles(gold_mol)
    # new_smiles = Chem.MolToSmiles(new_mol)
    # print("gold")
    # print(gold_smiles)
    # print("new")
    # print(new_smiles)
    # if gold_smiles != new_smiles:
    #     print('ERROR! The smiles are not the same.')
    return
        


def bond_to_label_2nd(bond):
    '''This function takes an RDKit bond and creates a label describing
    the most important attributes'''
    
    a1_label = str(bond.GetBeginAtom().GetAtomicNum())
    a2_label = str(bond.GetEndAtom().GetAtomicNum())
    if bond.GetBeginAtom().GetAtomMapNum():
        a1_label += str(bond.GetBeginAtom().GetIdx())
    if bond.GetEndAtom().GetAtomMapNum():
        a2_label += str(bond.GetEndAtom().GetIdx())
    atoms = sorted([a1_label, a2_label])

    return '{}{}'.format(atoms[0], atoms[1])


def atoms_are_different_2nd(atom1, atom2):
    '''Compares two RDKit atoms based on basic properties'''

    #if atom1.GetSmarts() != atom2.GetSmarts(): return True # should be very general
    if atom1.GetAtomicNum() != atom2.GetAtomicNum(): return True # must be true for atom mapping
    if atom1.GetTotalNumHs() != atom2.GetTotalNumHs(): return True
    if atom1.GetFormalCharge() != atom2.GetFormalCharge(): return True
    if atom1.GetDegree() != atom2.GetDegree(): return True
    if atom1.GetNumRadicalElectrons() != atom2.GetNumRadicalElectrons(): return True
    if atom1.GetIsAromatic() != atom2.GetIsAromatic(): return True 

    # Check bonds and nearest neighbor identity
    bonds1 = sorted([bond_to_label_2nd(bond) for bond in atom1.GetBonds()]) 
    bonds2 = sorted([bond_to_label_2nd(bond) for bond in atom2.GetBonds()]) 
    if bonds1 != bonds2: return True

    return False

if __name__ == '__main__':
    pass
    # example_reaction = '[CH3:1]/[CH:2]=[CH:3]/[CH:4]=[CH:5]/[CH2:6][O:7][C:8](=[O:9])[C:10]#[C:11][c:12]1[cH:13][cH:14][cH:15][cH:16][cH:17]1>>[CH3:1][C@@H:2]1[CH:3]=[CH:4][C@@H:5]2[CH2:6][O:7][C:8](=[O:9])[C:10]2=[C:11]1[c:12]1[cH:13][cH:14][cH:15][cH:16][cH:17]1'
    # reactants, product = example_reaction.split('>>')
    # reactants = ReassignMapping(reactants,iso=True)
    # product = ReassignMapping(product,iso=True)
    
    # smiles_list = rdchiralRunText_modified('[*:3]1-[*:4]~[*:5]-[*:6]-[*:1]~[*:2]-1>>[*:3]=[*:4]-[*:5]=[*:6].[*:1]#[*:2]',
    #                                        product)
    
    # new_reactants = smiles_list[0]
    # gold_mol = Chem.MolFromSmiles(reactants)
    # new_mol = Chem.MolFromSmiles(new_reactants)