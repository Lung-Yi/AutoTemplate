# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 11:40:40 2021

@author: Lung-Yi
"""

from rdkit import Chem
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(4)


BOND_TYPE = [0, Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC] 
BOND_FLOAT_TO_TYPE = {
    0.0: BOND_TYPE[0],
    1.0: BOND_TYPE[1],
    2.0: BOND_TYPE[2],
    3.0: BOND_TYPE[3],
    1.5: BOND_TYPE[4],
}
atom_valence_dict = {"C":4, "N":3, "O":2, "S":(2,4,6), "P":(3,5),"Si":4, "As":(3,5), "F":1, "Cl":1, "Br":1, "I":1,"B":3 }

def canon_remap(smiles):
    # copy mol before changing it
    mol = Chem.MolFromSmiles(smiles)
    if mol == None: return None
    # remove all atom numbers beforehand, as the CanonicalRankAtoms function
    # takes into account existing mapping numbers while doing it's mapping (!)
    for atom in mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            atom.ClearProp('molAtomMapNumber')

    return Chem.MolToSmiles(mol,isomericSmiles=False)

def CalculateNumHs(atom, valence):
    return max(0, int(valence - sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])) + atom.GetFormalCharge())

def UpdateAtomInfo(atom):
    """
    Checking the connectivity of the atom, if the explicit valence is beyond the maximum of atom valence,
    we set formal charge for the atom and then update the number of explicit hydrogens.
    """
    valence = atom_valence_dict.get(atom.GetSymbol())
    bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
    if type(valence) == int: 
        v_number = valence
    elif type(valence) == tuple:
        v_number = max(valence)
    if (atom.GetSymbol() == 'N') and (v_number < bond_vals):
        atom.SetFormalCharge(int(bond_vals) - v_number )
    
    if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 1: # exclude negatively-charged azide
        bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
        if bond_vals <= 3:
            atom.SetFormalCharge(0)
    elif atom.GetSymbol() == 'N' and atom.GetFormalCharge() == -1: # handle negatively-charged azide addition
        bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
        if bond_vals == 3 and any([nbr.GetSymbol() == 'N' for nbr in atom.GetNeighbors()]):
            atom.SetFormalCharge(0)
    elif atom.GetSymbol() == 'N':
        bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
        if bond_vals == 4 and not atom.GetIsAromatic(): # and atom.IsInRingSize(5)):
            atom.SetFormalCharge(1)
    elif atom.GetSymbol() == 'C' and atom.GetFormalCharge() != 0:
        atom.SetFormalCharge(0)
    elif atom.GetSymbol() == 'O' and atom.GetFormalCharge() != 0:
        bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])# + atom.GetNumExplicitHs()
        if bond_vals < 2: #(bond_vals == 2) or (bond_vals == 0):
            atom.SetFormalCharge(0)
    elif atom.GetSymbol() in ['Cl', 'Br', 'I', 'F'] and atom.GetFormalCharge() != 0:
        bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
        if bond_vals == 1:
            atom.SetFormalCharge(0)
    elif atom.GetSymbol() == 'S' and atom.GetFormalCharge() != 0:
        bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
        if bond_vals in [2, 4, 6]:
            atom.SetFormalCharge(0)
    elif atom.GetSymbol() == 'P': # quartenary phosphorous should be pos. charge with 0 H
        bond_vals = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
        if sum(bond_vals) == 4 and len(bond_vals) == 4:
            atom.SetFormalCharge(1)
            atom.SetNumExplicitHs(0)
        elif sum(bond_vals) == 3 and len(bond_vals) == 3: # make sure neutral
            atom.SetFormalCharge(0)
    elif atom.GetSymbol() == 'B': # quartenary boron should be neg. charge with 0 H
        bond_vals = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
        if sum(bond_vals) == 4 and len(bond_vals) == 4:
            atom.SetFormalCharge(-1)
            atom.SetNumExplicitHs(0)
    elif atom.GetSymbol() in ['Mg', 'Zn']:
        bond_vals = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
        if sum(bond_vals) == 1 and len(bond_vals) == 1:
            atom.SetFormalCharge(1)
    elif atom.GetSymbol() == 'Si':
        bond_vals = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
        if sum(bond_vals) == len(bond_vals):
            atom.SetNumExplicitHs(max(0, 4 - len(bond_vals)))
    
    if type(valence) == int:
        try:
            atom.SetNumExplicitHs(CalculateNumHs(atom, valence))
            atom.UpdatePropertyCache(strict = True) # test: strict
        except:
            return
            # print(Chem.MolToSmiles(atom.GetOwningMol()))
            # print(atom.GetSmarts())
    elif type(valence) == tuple:
        for v in valence:
            try:
                hydrogens = CalculateNumHs(atom, v)
                if hydrogens < 0: continue
                atom.SetNumExplicitHs(hydrogens)
                atom.UpdatePropertyCache(strict = True) # test: strict
                break
            except:
                print(Chem.MolToSmiles(atom.GetOwningMol()))
                print(atom.GetSmarts())
                continue
    return

def edit_mol(rmol, edits): # 2021/02/25: change inside function
    """
    return chiral, cis-trans informaton and atom-mapping product
    """
    new_mol = Chem.RWMol(rmol)
    Chem.Kekulize(new_mol)
    amap = {}
    for atom in rmol.GetAtoms():
        amap[atom.GetAtomMapNum()] = atom.GetIdx()
    
    changed_atom_idx = set()
    for x,y,t in edits:
        changed_atom_idx.add(amap[x])
        changed_atom_idx.add(amap[y])
        bond = new_mol.GetBondBetweenAtoms(amap[x],amap[y])
        # a1 = new_mol.GetAtomWithIdx(amap[x])
        # a2 = new_mol.GetAtomWithIdx(amap[y])
        if bond is not None:
            new_mol.RemoveBond(amap[x],amap[y])
        if t > 0:
            new_mol.AddBond(amap[x],amap[y],BOND_FLOAT_TO_TYPE[t])
    for idx in changed_atom_idx:
        change_atom = new_mol.GetAtomWithIdx(idx)
        UpdateAtomInfo(change_atom)
        change_atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
        change_atom.ClearProp('_CIPCode')
        # z = -3 # change formal charge from -3 to 4
        # while True:
        #     try:
        #         change_atom.UpdatePropertyCache(strict=True)
        #         break
        #     except Chem.AtomValenceException:
        #         change_atom.SetFormalCharge(z)
        #         z += 1
        #         if z ==5 : break
    for idx in changed_atom_idx:
        change_atom = new_mol.GetAtomWithIdx(idx)
        z = -3 # change formal charge from -3 to 4
        while True:
            try:
                change_atom.UpdatePropertyCache(strict=True)
                break
            except Chem.AtomValenceException:
                change_atom.SetFormalCharge(z)
                z += 1
                if z ==3 : break
    new_mol.UpdatePropertyCache(strict=True)
    try:
        Chem.SanitizeMol(new_mol)
        Chem.SetAromaticity(new_mol) # set Aromaticity 
    except:
        for atom in new_mol.GetAtoms():
            if not atom.IsInRing():
                atom.SetIsAromatic(False)
        try:
            Chem.SanitizeMol(new_mol)
            Chem.SetAromaticity(new_mol)
        except:
            return
    
    chiral_index = Chem.FindMolChiralCenters(new_mol, includeUnassigned=True)
    for idx, tag in chiral_index:
        if tag == '?':
            pass
            # print('potential chiral atom: {}'.format(new_mol.GetAtomWithIdx(idx).GetAtomMapNum()))
    # pred_mol = new_mol.GetMol()
    pred_smiles = Chem.MolToSmiles(new_mol)
    return pred_smiles