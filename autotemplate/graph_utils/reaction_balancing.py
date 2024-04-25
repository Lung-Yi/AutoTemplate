import CGRtools
from CGRtools.containers import ReactionContainer, MoleculeContainer
from rdkit import Chem
import numpy as np
import sys
sys.path.append("../../")
from autotemplate.run_utils import RemoveReagent

def give_mapnum_for_unmapped_reactant_atoms(rxn_smiles):
    """ 
    Leaving group atoms in reactnat site have no atom mapping number.
    Try to append atom mapping number for those reactant atoms.
    """
    reac, prod = rxn_smiles.split('>>')
    reac = Chem.MolFromSmiles(reac)
    prod = Chem.MolFromSmiles(prod)
    if reac.GetNumAtoms() > prod.GetNumAtoms():
        reac_map = {i+1 for i in range(reac.GetNumAtoms())}
        prod_map = {int(atom.GetProp('molAtomMapNumber')) for atom in prod.GetAtoms()}
        undetermined_map = reac_map.difference(prod_map)
        for atom in reac.GetAtoms():
            if not atom.HasProp('molAtomMapNumber'):
                num = undetermined_map.pop()
                atom.SetAtomMapNum(num)
                if not undetermined_map: break
        return Chem.MolToSmiles(reac) + '>>' + Chem.MolToSmiles(prod)
    else:
        return rxn_smiles

def SmilesFromGraphs(node_list, adjacency_matrix, mapping_num_list, charge_list):
    mol = Chem.RWMol()
    node_to_idx = {}
    idx_to_mapping_num = {}
    idx_to_charge = {}
    for i in range(len(node_list)):
        a = Chem.Atom(node_list[i])
        molIdx = mol.AddAtom(a)
        node_to_idx[i] = molIdx
        idx_to_mapping_num[molIdx] = mapping_num_list[i]
        idx_to_charge[molIdx] = charge_list[i]
    # add bonds between adjacent atoms
    for ix, row in enumerate(adjacency_matrix):
        for iy, bond in enumerate(row):
            # only traverse half the matrix
            if iy <= ix:
                continue
            if bond == 0:
                continue
            elif bond == 1:
                bond_type = Chem.rdchem.BondType.SINGLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
            elif bond == 2:
                bond_type = Chem.rdchem.BondType.DOUBLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
            elif bond == 3:
                bond_type = Chem.rdchem.BondType.TRIPLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
            elif bond == 4:
                bond_type = Chem.rdchem.BondType.AROMATIC
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
    mol = mol.GetMol()
    for molIdx, mapping_num in idx_to_mapping_num.items():
        atom = mol.GetAtomWithIdx(molIdx)
        atom.SetAtomMapNum(mapping_num)
        atom.SetFormalCharge(idx_to_charge[molIdx])
    return Chem.MolToSmiles(mol)

def ConnectivityFromContainer(mc: MoleculeContainer):
    node_list = []
    mapping_num_list = []
    charge_list = []
    mapping_num_to_idx = {}
    for i, atom in enumerate(mc.atoms()):
        map_num, atom_info = atom
        node_list.append(atom_info.atomic_symbol)
        charge_list.append(atom_info.charge)
        mapping_num_list.append(map_num)
        mapping_num_to_idx.update({map_num: i})

    adjacency_matrix = np.zeros((len(node_list), len(node_list)))
    for bond in mc.bonds():
        map_num_1, map_num_2, bond_info = bond
        bond_order = bond_info.order
        i = mapping_num_to_idx[map_num_1]
        j = mapping_num_to_idx[map_num_2]
        adjacency_matrix[i][j] = bond_order
        adjacency_matrix[j][i] = bond_order
    return node_list, adjacency_matrix, mapping_num_list, charge_list

def CleanStereo(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles),isomericSmiles=False)

def balance_rxn_smiles(rxn_smiles):
    """ This function only balances the product site, since the reactant site has been balanced. """
    rxn_smiles = give_mapnum_for_unmapped_reactant_atoms(rxn_smiles)
    cleaned_rxn_smiles, reagent_smiles = RemoveReagent(rxn_smiles, retain_reagents=True)
    reactant, product = cleaned_rxn_smiles.split(">>")
    if Chem.MolFromSmiles(reactant).GetNumAtoms() <= Chem.MolFromSmiles(product).GetNumAtoms():
        return rxn_smiles
    r = CGRtools.smiles(reactant)
    p = CGRtools.smiles(product)
    cgr = r^p

    decomposed = ReactionContainer.from_cgr(cgr)
    decomposed.clean2d()
    output_rxn_smiles = ""
    for reac_mc in decomposed.reactants:
        node_list, adjacency_matrix, mapping_num_list, charge_list = ConnectivityFromContainer(reac_mc)
        smiles = SmilesFromGraphs(node_list, adjacency_matrix, mapping_num_list, charge_list)
        output_rxn_smiles = output_rxn_smiles + smiles + "."
    output_rxn_smiles = output_rxn_smiles.strip(".")
    output_rxn_smiles = output_rxn_smiles + ">>"
    
    original_p = CleanStereo(p.__str__())
    new_p_mcs = [x for x in decomposed.products if CleanStereo(x.__str__()) != original_p]
    output_rxn_smiles = cleaned_rxn_smiles
    if new_p_mcs == []: # no new molecules generated after reaction balancing
        return rxn_smiles
    for p_mc in new_p_mcs:
        node_list, adjacency_matrix, mapping_num_list, charge_list = ConnectivityFromContainer(p_mc)
        smiles = SmilesFromGraphs(node_list, adjacency_matrix, mapping_num_list, charge_list)
        output_rxn_smiles = output_rxn_smiles + "." + smiles
    if reagent_smiles:
        output_rxn_smiles = reagent_smiles + "." + output_rxn_smiles
    return output_rxn_smiles


if __name__ == "__main__":
    # Reaction taken from the USPTO-50k dataset, with reaction ID: US20100009970A1_44
    rxn_smiles = "CC(C)(C)OC(=O)O[C:6]([O:5][C:2]([CH3:1])([CH3:3])[CH3:4])=[O:7].[CH3:8][NH:9][C@H:10]1[CH2:11][CH2:12][C@@H:13]([c:14]2[cH:15][cH:16][c:17]([Cl:18])[c:19]([Cl:20])[cH:21]2)[c:22]2[cH:23][cH:24][c:25]([C:26](=[O:27])[O:28][CH3:29])[cH:30][c:31]21>>[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:9]([CH3:8])[C@H:10]1[CH2:11][CH2:12][C@@H:13]([c:14]2[cH:15][cH:16][c:17]([Cl:18])[c:19]([Cl:20])[cH:21]2)[c:22]2[cH:23][cH:24][c:25]([C:26](=[O:27])[O:28][CH3:29])[cH:30][c:31]21"
    print(balance_rxn_smiles(rxn_smiles))