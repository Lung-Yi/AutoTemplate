from rdkit import Chem
from rdchiral.template_extractor import mols_from_smiles_list

def get_reaction_atom_difference(rxn_smiles):
    """ Calculate the number difference of reactant and product site. """
    reac, prod = rxn_smiles.split('>>')
    reac_list = reac.split('.')
    prod_list = prod.split('.')
    reac_mol_list = mols_from_smiles_list(reac_list)
    prod_mol_list = mols_from_smiles_list(prod_list)
    return sum([mol.GetNumAtoms() for mol in prod_mol_list]) - sum([mol.GetNumAtoms() for mol in reac_mol_list])

if __name__ == "__main__":
    rxn_smiles = "[CH2:1]=[CH:10][CH:11]=[O:12].[CH:2]1=[CH:3][CH:9]=[CH:5][CH2:6][CH2:7][CH2:8]1>>[CH3:1][C@:2]12[CH:3]=[CH:4][C@H:5]([CH2:6][CH2:7][CH2:8]1)[CH2:9][C@H:10]2[CH:11]=[O:12]"
    print(get_reaction_atom_difference(rxn_smiles))