from rdkit import Chem
from rdkit.Chem import AllChem
from rdchiral.main import rdchiralRunText
from rdchiral.utils import PLEVEL
from rdchiral.initialization import rdchiralReaction, rdchiralReactants

atom_valence_dict = {"C":4, "N":3, "O":2, "S":(2,4,6), "P":(3,5),"Si":4, "Br":1, "Cl":1, "I":1, "F":1, "Mg":2}

def clearIsotope(rxn_smiles):
    reac, prod = rxn_smiles.split('>>')
    m_r = Chem.MolFromSmiles(reac)
    m_p = Chem.MolFromSmiles(prod)
    [atom.SetIsotope(0) for atom in m_r.GetAtoms()]
    [atom.SetIsotope(0) for atom in m_p.GetAtoms()]
    return Chem.MolToSmiles(m_r) + '>>' + Chem.MolToSmiles(m_p)

def CalculateNumHs(atom):
    valence = atom_valence_dict[atom.GetSymbol()]
    return int(valence - sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]) - abs(atom.GetFormalCharge()))

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

if __name__ == "__main__":
    test_reaction = "[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[NH:8][OH:9].[CH:10]1=[CH:11][CH:12]=[CH:13][CH2:14][CH2:15]1>>[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[O:9][C@H:10]2[CH:11]=[CH:12][C@@H:13]1[CH2:14][CH2:15]2"
    template_0 = "[#6:1]1-[#6:2]=[#6:3]-[#6:4]-[#7:5]-[#8:6]-1>>[#6:1]=[#6:2]-[#6:3]=[#6:4].[#7:5]-[#8:6]"
    template_1 = "[#6:9]-[#7:8]1-[#8:7]-[#6:1]2-[#6:2]=[#6:3]-[#6:4]-1-[#6:5]-[#6:6]-2>>[#6:1]1=[#6:2]-[#6:3]=[#6:4]-[C:5]-[C:6]-1.[#8:7]-[#7:8]-[C:9]"
    reactants, product = test_reaction.split(">>")
    print("test_reaction")
    print("template radius 0:")
    print(template_0)
    print(rdchiralRunText_modified(template_0, product))
    print("template radius 1:")
    print(template_1)
    print(rdchiralRunText_modified(template_1, product))
    # rdchiralRunText(reaction_smarts, reactant_smiles, keep_mapnums=True, return_mapped=True)
    