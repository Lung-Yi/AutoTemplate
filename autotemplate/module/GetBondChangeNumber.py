# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:45:56 2020

@author: Lung-Yi
"""
from rdkit import Chem
# from .edit_mol_function import edit_mol, canon_remap

#rxn_smiles = '[CH2:15]([CH:16]([CH3:17])[CH3:18])[Mg+:19].[CH2:20]1[O:21][CH2:22][CH2:23][CH2:24]1.[Cl-:14].[OH:1][c:2]1[n:3][cH:4][c:5]([C:6](=[O:7])[N:8]([O:9][CH3:10])[CH3:11])[cH:12][cH:13]1>>[OH:1][c:2]1[n:3][cH:4][c:5]([C:6](=[O:7])[CH2:15][CH:16]([CH3:17])[CH3:18])[cH:12][cH:13]1'
#answer = '6-8-0.0;15-6-1.0;15-19-0.0'

def swap(i,j):
    return min(i,j), max(i,j)

def DeleteTrivialProduct(rxn_smiles):
    '''Delete trivial molecules that have number of heavy atoms less than 5'''
    rea = rxn_smiles.split('>')[0]
    prod = rxn_smiles.split('>')[-1]
    if '.' in prod:
        prod_list = prod.split('.')
        prod_list2 = prod_list[:]
        for p in prod_list:
            if Chem.MolFromSmiles(p).GetNumAtoms() < 5:
                prod_list2.remove(p)
        prod = '.'.join(prod_list2)
        return rea + '>>' + prod
    else:
        return rxn_smiles

def CheckMassBalance(rxn_smiles):
    rea = rxn_smiles.split('>')[0]
    prod = rxn_smiles.split('>')[-1]
    r = Chem.MolFromSmiles(rea)
    p = Chem.MolFromSmiles(prod)
    r_info = {}
    p_info = {}
    for atom in r.GetAtoms():
        r_info.update({atom.GetAtomMapNum():atom})
    for atom in p.GetAtoms():
        p_info.update({atom.GetAtomMapNum():atom})
    set_r = set(r_info.keys())
    set_p = set(p_info.keys())
    return set_p.issubset(set_r)


def AppendAtomMapNum(rxn_smiles):
    rea = rxn_smiles.split('>')[0]
    mol = Chem.MolFromSmiles(rea)
    atom_index = []
    for atom in mol.GetAtoms():
        num = atom.GetAtomMapNum()
        if num != 0:
            atom_index.append(num)
    
    i = mol.GetNumAtoms()
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            atom.SetAtomMapNum(i)
            i -= 1
            while i in atom_index:
                i -= 1
    rea = Chem.MolToSmiles(mol)
    return rea + '>>' + rxn_smiles.split('>>')[1]

def CheckZeroIndex(rxn_smiles):
    '''cannot have zero index atom-mapping number in smiles, or the model cannot be trained'''
    rea = rxn_smiles.split('>>')[0]
    prod = rxn_smiles.split('>>')[-1]
    rea_mol = Chem.MolFromSmiles(rea)
    prod_mol = Chem.MolFromSmiles(prod)
    for atom in rea_mol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            return True
    for atom in prod_mol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            return True
    return False

def CheckDuplicateIndex_andContinuity(rxn_smiles):
    '''cannot have zero index atom-mapping number in smiles, or the model cannot be trained
    Check if all the atom-mapping index is continous number'''
    rea = rxn_smiles.split('>>')[0]
    prod = rxn_smiles.split('>>')[1]
    rea_mol = Chem.MolFromSmiles(rea)
    prod_mol = Chem.MolFromSmiles(prod)
    if rea_mol == None:
        return True
    if prod_mol == None:
        return True
    rea_index = []
    for atom in rea_mol.GetAtoms():
        num = atom.GetAtomMapNum()
        if num in rea_index:
            return True
        rea_index.append(num)
    return not (CheckContinuity(rea_index))

    return False

def CheckContinuity(index:list):
    ans = [x+1 for x in range(len(index))]
    for x in index:
        if x not in ans:
            return False
    return True


def GetAnswer(rxn_smiles, kekulize = True):
    AllBondChange = []
    rea = rxn_smiles.split('>')[0]
    prod = rxn_smiles.split('>')[-1]
    r = Chem.MolFromSmiles(rea)
    p = Chem.MolFromSmiles(prod)
    if kekulize:
        Chem.Kekulize(r) # must have this to solve aromatic ring problem, such as pyrrole
        Chem.Kekulize(p) #
    
    r_info = {}
    p_info = {}
    for atom in r.GetAtoms():
        r_info.update({atom.GetAtomMapNum():atom})
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


def GetMixAnswer_New(rxn_smiles):
    """
    Better Get Bond Change Answer function
    """
    AllBondChange = []
    rea = rxn_smiles.split('>')[0]
    prod = rxn_smiles.split('>')[-1]
    r_sani = Chem.MolFromSmiles(rea)
    p_sani = Chem.MolFromSmiles(prod)
    
    r_kekule = r_sani.__copy__()
    p_kekule = p_sani.__copy__()
    Chem.Kekulize(r_kekule) # must have this to solve aromatic ring problem, such as pyrrole
    Chem.Kekulize(p_kekule) #
    
    r_sani_info = {}
    p_sani_info = {}
    r_kekule_info = {}
    p_kekule_info = {}
    for atom in r_sani.GetAtoms():
        r_sani_info.update({atom.GetAtomMapNum():atom.GetIdx()})
    for atom in p_sani.GetAtoms():
        p_sani_info.update({atom.GetAtomMapNum():atom.GetIdx()})
    for atom in r_kekule.GetAtoms():
        r_kekule_info.update({atom.GetAtomMapNum():atom.GetIdx()})
    for atom in p_kekule.GetAtoms():
        p_kekule_info.update({atom.GetAtomMapNum():atom.GetIdx()})
    
    for mapnum, idx in p_kekule_info.items(): # check first kekule form, then aromatic form
        atom_r_kekule = r_kekule.GetAtomWithIdx( r_kekule_info.get(mapnum) )
        atom_p_kekule = p_kekule.GetAtomWithIdx( idx )
        result_kekule = CheckTwoAtoms(atom_r_kekule, atom_p_kekule)
        
        if result_kekule:
            atom_r_sani = r_sani.GetAtomWithIdx( r_sani_info.get(mapnum) )
            atom_p_sani = p_sani.GetAtomWithIdx( p_sani_info.get(mapnum) )
            result_sani = CheckTwoAtoms(atom_r_sani, atom_p_sani)
            if result_sani:
                changes = ['-'.join(result.split('-')[:-1]) for result in result_sani]
                true_result = []
                for change in changes:
                    for result in result_kekule:
                        if change in result: 
                            true_result.append(result)
                AllBondChange += true_result
        
    AllBondChange = list(set(AllBondChange))
    return ';'.join(AllBondChange)


def GetMixAnswer_New_Beta(rxn_smiles):
    """
    Get Bond Change Answer function. Considering the enumeration of different kekule form of product.
    """
    
    rea = rxn_smiles.split('>')[0]
    prod = rxn_smiles.split('>')[-1]
    r_sani = Chem.MolFromSmiles(rea)
    p_sani = Chem.MolFromSmiles(prod)
    
    r_kekule = r_sani.__copy__()
    p_kekule_raw = p_sani.__copy__()
    Chem.Kekulize(r_kekule) # must have this to solve aromatic ring problem, such as pyrrole
    Chem.Kekulize(p_kekule_raw) #
    p_kekule_list = enumerate_kekule_mol(p_kekule_raw)
    
    BCN = 20 # maximum Bond change number
    DEBUG = []
    for p_kekule in p_kekule_list:
        AllBondChange = []
        r_sani_info = {}
        p_sani_info = {}
        r_kekule_info = {}
        p_kekule_info = {}
        for atom in r_sani.GetAtoms():
            r_sani_info.update({atom.GetAtomMapNum():atom.GetIdx()})
        for atom in p_sani.GetAtoms():
            p_sani_info.update({atom.GetAtomMapNum():atom.GetIdx()})
        for atom in r_kekule.GetAtoms():
            r_kekule_info.update({atom.GetAtomMapNum():atom.GetIdx()})
        for atom in p_kekule.GetAtoms():
            p_kekule_info.update({atom.GetAtomMapNum():atom.GetIdx()})
        
        for mapnum, idx in p_kekule_info.items(): # check first kekule form, then aromatic form
            atom_r_kekule = r_kekule.GetAtomWithIdx( r_kekule_info.get(mapnum) )
            atom_p_kekule = p_kekule.GetAtomWithIdx( idx )
            result_kekule = CheckTwoAtoms(atom_r_kekule, atom_p_kekule)
            
            if result_kekule:
                # atom_r_sani = r_sani.GetAtomWithIdx( r_sani_info.get(mapnum) )
                # atom_p_sani = p_sani.GetAtomWithIdx( p_sani_info.get(mapnum) )
                # result_sani = CheckTwoAtoms(atom_r_sani, atom_p_sani)
                # if result_sani:
                #     changes = ['-'.join(result.split('-')[:-1]) for result in result_sani]
                #     true_result = []
                #     for change in changes:
                #         for result in result_kekule:
                #             if change in result: 
                #                 true_result.append(result)
                AllBondChange += result_kekule# true_result
        
        AllBondChange = list(set(AllBondChange))
        DEBUG.append(AllBondChange)
        if len(AllBondChange) < BCN:
            if IsAnswerValid(r_kekule, prod, AllBondChange):
                BCN = len(AllBondChange)
                Kekule_answer = ';'.join(AllBondChange)
            
    if BCN == 20:
        return "Answer Not Found"
    else:
        return Kekule_answer

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
    
def CountNumBondChange(answer):
    if ';' not in answer:
        return 1
    else:
        return answer.count(';')+1

def IsAnswerValid(rmol, psmi, AllBondChange):
    """
    rmol: Chem.Mol needed to be kekule form
    psmi: Chem.Mol Aromatic smiles
    Check the answer calculated by GetAnswer function is appliable for edit_mol function.
    If the product can not be matched with the sudo-product, then the answer if False.
    """
    gold_edits = []
    for edit in AllBondChange:
        x, y, t = edit.split('-')
        g = (int(x), int(y), float(t))
        gold_edits.append(g)
        
    try_psmi = edit_mol(rmol, gold_edits)
    if try_psmi == None: 
        return False
    try_psmi_copy = canon_remap(try_psmi)
    if try_psmi_copy == None:
        return False
    
    try_psmi_copy = try_psmi_copy.split('.')
    if [i for i in canon_remap(psmi).split('.') if not i in try_psmi_copy or try_psmi_copy.remove(i)]:
        return False
    return True
    


"""
Enumerate ALL possible Kekule form of molecule
"""
MAX_STRUCTURES = 1000

class ResonanceEnumerator(object):
    """Simple wrapper around RDKit ResonanceMolSupplier.
    """

    def __init__(self, kekule_all=False, allow_incomplete_octets=False, unconstrained_cations=True,
                 unconstrained_anions=False, allow_charge_separation=True, max_structures=MAX_STRUCTURES):
        """
        :param bool allow_incomplete_octets: include resonance structures whose octets are less complete than the the most octet-complete structure.
        :param bool allow_charge_separation: include resonance structures featuring charge separation also when uncharged resonance structures exist.
        :param bool kekule_all: enumerate all possible degenerate Kekule resonance structures (the default is to include just one).
        :param bool unconstrained_cations: if False positively charged atoms left and right of N with an incomplete octet are acceptable only if the conjugated group has a positive total formal charge.
        :param bool unconstrained_anions: if False, negatively charged atoms left of N are acceptable only if the conjugated group has a negative total formal charge.
        :param int max_structures: Maximum number of resonance forms.
        """
        self.kekule_all = kekule_all
        self.allow_incomplete_octets = allow_incomplete_octets
        self.unconstrained_cations = unconstrained_cations
        self.unconstrained_anions = unconstrained_anions
        self.allow_charge_separation = allow_charge_separation
        self.max_structures = max_structures

    def __call__(self, mol):
        """Calling a ResonanceEnumerator instance like a function is the same as calling its enumerate(mol) method."""
        return self.enumerate(mol)

    def enumerate(self, mol):
        """Enumerate all possible resonance forms and return them as a list.
        :param mol: The input molecule.
        :type mol: rdkit.Chem.rdchem.Mol
        :return: A list of all possible resonance forms of the molecule.
        :rtype: list of rdkit.Chem.rdchem.Mol
        """
        flags = 0
        if self.kekule_all:
            flags = flags | Chem.KEKULE_ALL
        if self.allow_incomplete_octets:
            flags = flags | Chem.ALLOW_INCOMPLETE_OCTETS
        if self.allow_charge_separation:
            flags = flags | Chem.ALLOW_CHARGE_SEPARATION
        if self.unconstrained_anions:
            flags = flags | Chem.UNCONSTRAINED_ANIONS
        if self.unconstrained_cations:
            flags = flags | Chem.UNCONSTRAINED_CATIONS
        results = []
        for result in Chem.ResonanceMolSupplier(mol, flags=flags, maxStructs=self.max_structures):
            # This seems necessary? ResonanceMolSupplier only does a partial sanitization
            # Chem.SanitizeMol(result)
            results.append(result)
        return results
    
def enumerate_kekule_mol(mol):
    """Return a set of resonance forms as SMILES strings, given a SMILES string.
    :param smiles: A SMILES string.
    :returns: A set containing SMILES strings for every possible resonance form.
    :rtype: set of strings.
    """
    Chem.SanitizeMol(mol)
    Chem.Kekulize(mol)
    mesomers = ResonanceEnumerator(kekule_all=True).enumerate(mol)
    return mesomers

"""
Data preprocess for USPTO-MIT dataset
"""

def AtomMap_reagent(smiles):
    mol = Chem.MolFromSmiles(smiles)
    number = mol.GetNumAtoms()
    for atom in mol.GetAtoms():
        if not atom.HasProp('molAtomMapNumber'):
            atom.SetAtomMapNum( number)
            number -= 1
    return Chem.MolToSmiles(mol, isomericSmiles=True)



# if __name__ == "__main__":
#     from tqdm import tqdm
#     import os
#     data_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#     data_path = os.path.join(os.path.join(data_path, 'data'), '11template_splitted')
#     data_path = os.path.join(os.path.join(data_path, 'With_StereoChiral_mix'), 'filtered_Hydrogenation.txt')
#     # f = open('D:\\Retro\\LCC_project\\data\\10template_processed\\With_StereoChiral_mix\\TPed_Diels_Alder_ChiralStereo.txt')
#     f = open(data_path, 'r')
#     data = f.readlines()
#     f.close()
#     out_path = os.path.join(os.path.dirname(os.path.dirname(data_path)), 'filtered_Hydrogenation_ZZ.txt')
#     g = open(out_path, 'w')
    
#     count_list = []
#     for line in tqdm(data):
#         rxn_smiles, answer_ = line.strip('\n').split(' ')
#         answer = GetMixAnswer_New_Beta(rxn_smiles)
#         CCC = CountNumBondChange(answer)
#         count_list.append(CCC)
#         if (CCC != 1) and (CCC<=6):
#             g.write(rxn_smiles+' '+answer+'\n')
        
#     from collections import Counter
#     results = Counter(count_list)
#     g.close()
#     print(results)
