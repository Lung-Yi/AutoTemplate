from rdchiral.template_extractor import mols_from_smiles_list, replace_deuterated, get_changed_atoms, clear_mapnum, get_tagged_atoms_from_mols, \
                            get_tagged_atoms_from_mol, get_fragments_for_changed_atoms, \
                            find_map_num, get_tetrahedral_atoms, set_isotope_to_equal_mapnum, get_frag_around_tetrahedral_center, check_tetrahedral_centers_equivalent, \
                            clear_isotope, expand_atoms_to_use, expand_atoms_to_use_atom, reassign_atom_mapping, expand_changed_atom_tags, canonicalize_transform, \
                            canonicalize_template

from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.rdchem import ChiralType
import re
from numpy.random import shuffle
from copy import deepcopy
import networkx as nx

VERBOSE = False
USE_STEREOCHEMISTRY = True
MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS = 5
INCLUDE_ALL_UNMAPPED_REACTANT_ATOMS = True

def canon_remap(smiles, return_NumAtom=False, iso = False):
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
    return Chem.MolToSmiles(mol,isomericSmiles=iso)

def get_atoms_attributes_template(rdkit_molecule):
    attributes = []
    for a in rdkit_molecule.GetAtoms():
        attributes.append((a.GetAtomMapNum(), {'Idx':a.GetIdx()}))
    return attributes

def get_bonds_attributes_template(rdkit_molecule):
    attributes = []
    for b in rdkit_molecule.GetBonds():
        start = b.GetBeginAtom().GetAtomMapNum()
        end = b.GetEndAtom().GetAtomMapNum()
        attributes.append((start, end))
    return attributes

def topology_template_from_rdkit(rdkit_molecule):
    topology = nx.Graph()
    topology.add_nodes_from(get_atoms_attributes_template(rdkit_molecule))
    topology.add_edges_from(get_bonds_attributes_template(rdkit_molecule))
    return topology

def HaveSeparateReactionCenter(rdkit_molecule, changed_atom_tags):
    if not changed_atom_tags: return False
    G = topology_template_from_rdkit(rdkit_molecule)
    changed_atom_tags = [int(x) for x in changed_atom_tags]
    # print(Chem.MolToSmiles(rdkit_molecule))
    # print(changed_atom_tags)
    H = G.subgraph(changed_atom_tags)
    return not nx.is_connected(H)

def UpdateChangedAtomTags(rdkit_molecule, changed_atom_tags):
    G = topology_template_from_rdkit(rdkit_molecule)
    changed_atom_tags = [int(x) for x in changed_atom_tags]
    H = G.subgraph(changed_atom_tags)
    reaction_center_clusters = [c for c in nx.connected_components(H)]
    path_atoms = []
    for i in range(len(reaction_center_clusters) - 1):
        for j in range(i+1, len(reaction_center_clusters)):
            c1 = next(iter(reaction_center_clusters[i]))
            c2 = next(iter(reaction_center_clusters[j]))
            try:# The reaction centers maybe lie in different product molecules
                path = nx.dijkstra_path(G,c1,c2)
            except:
                continue
            for x in path:
                if x not in changed_atom_tags: 
                    changed_atom_tags.append(x)
                    path_atoms.append(G.nodes[x]['Idx'])
    
    changed_atom_tags = sorted(changed_atom_tags)
    return [str(x) for x in changed_atom_tags], path_atoms

def get_changed_atoms(reactants, products):
    '''Looks at mapped atoms in a reaction and determines which ones changed'''

    err = 0
    prod_atoms, prod_atom_tags = get_tagged_atoms_from_mols(products)

    if VERBOSE: print('Products contain {} tagged atoms'.format(len(prod_atoms)))
    if VERBOSE: print('Products contain {} unique atom numbers'.format(len(set(prod_atom_tags))))

    reac_atoms, reac_atom_tags = get_tagged_atoms_from_mols(reactants)
    if len(set(prod_atom_tags)) != len(set(reac_atom_tags)):
        if VERBOSE: print('warning: different atom tags appear in reactants and products')
        #err = 1 # okay for Reaxys, since Reaxys creates mass
    if len(prod_atoms) != len(reac_atoms):
        if VERBOSE: print('warning: total number of tagged atoms differ, stoichometry != 1?')
        #err = 1

    # Find differences 
    changed_atoms = [] # actual reactant atom species
    changed_atom_tags = [] # atom map numbers of those atoms

    # Product atoms that are different from reactant atom equivalent
    for i, prod_tag in enumerate(prod_atom_tags):

        for j, reac_tag in enumerate(reac_atom_tags):
            if reac_tag != prod_tag: continue
            if reac_tag not in changed_atom_tags: # don't bother comparing if we know this atom changes
                # If atom changed, add
                if atoms_are_different(prod_atoms[i], reac_atoms[j]):
                    changed_atoms.append(reac_atoms[j])
                    changed_atom_tags.append(reac_tag)
                    break
                # If reac_tag appears multiple times, add (need for stoichometry > 1)
                if prod_atom_tags.count(reac_tag) > 1:
                    changed_atoms.append(reac_atoms[j])
                    changed_atom_tags.append(reac_tag)
                    break
        # TODO: include the atoms only appear in product site but not in reactant
        # else:
            

    # Reactant atoms that do not appear in product (tagged leaving groups)
    for j, reac_tag in enumerate(reac_atom_tags):
        if reac_tag not in changed_atom_tags:
            if reac_tag not in prod_atom_tags:
                changed_atoms.append(reac_atoms[j])
                changed_atom_tags.append(reac_tag)
    
    
    # Atoms that change CHIRALITY (just tetrahedral for now...)
    tetra_atoms = get_tetrahedral_atoms(reactants, products)
    if VERBOSE:
        print('Found {} atom-mapped tetrahedral atoms that have chirality specified at least partially'.format(len(tetra_atoms)))
    [set_isotope_to_equal_mapnum(reactant) for reactant in reactants]
    [set_isotope_to_equal_mapnum(product) for product in products]
    for (atom_tag, ar, ap) in tetra_atoms:
        if VERBOSE: 
            print('For atom tag {}'.format(atom_tag))
            print('    reactant: {}'.format(ar.GetChiralTag()))
            print('    product:  {}'.format(ap.GetChiralTag()))
        if atom_tag in changed_atom_tags:
            if VERBOSE:
                print('-> atoms have changed (by more than just chirality!)')
        else:
            unchanged = check_tetrahedral_centers_equivalent(ar, ap) and \
                    ChiralType.CHI_UNSPECIFIED not in [ar.GetChiralTag(), ap.GetChiralTag()]
            if unchanged:
                if VERBOSE: 
                    print('-> atoms confirmed to have same chirality, no change')
            else:
                if VERBOSE:
                    print('-> atom changed chirality!!')
                # Make sure chiral change is next to the reaction center and not
                # a random specifidation (must be CONNECTED to a changed atom)
                tetra_adj_to_rxn = False
                for neighbor in ap.GetNeighbors():
                    if neighbor.HasProp('molAtomMapNumber'):
                        if neighbor.GetProp('molAtomMapNumber') in changed_atom_tags:
                            tetra_adj_to_rxn = True
                            break
                if tetra_adj_to_rxn:
                    if VERBOSE:
                        print('-> atom adj to reaction center, now included')
                    changed_atom_tags.append(atom_tag)
                    changed_atoms.append(ar)
                else:
                    if VERBOSE:
                        print('-> adj far from reaction center, not including')
    [clear_isotope(reactant) for reactant in reactants]
    [clear_isotope(product) for product in products]

    # TODO: Include all atoms between all the reaction center:
    p_smiles = '.'.join([Chem.MolToSmiles(p) for p in products])
    prod_mol = Chem.MolFromSmiles(p_smiles)
    if HaveSeparateReactionCenter(prod_mol, changed_atom_tags):
        if VERBOSE: print('This template in product site has separate reaction centers.')
        # print(changed_atom_tags)
        changed_atom_tags, path_atoms = UpdateChangedAtomTags(prod_mol, changed_atom_tags)
        # print(changed_atom_tags)
  
    if VERBOSE: 
        print('{} tagged atoms in reactants change 1-atom properties'.format(len(changed_atom_tags)))
        for smarts in [atom.GetSmarts() for atom in changed_atoms]:
            print('  {}'.format(smarts))

    return changed_atoms, changed_atom_tags, err

def convert_atom_to_wildcard(atom):
    '''This function is modified. Only the most general information is kept. (Atom number)
    The degree of hydrogen is removed. This is generic reaction template version.'''
    
    if ':' in atom.GetSmarts():
        symbol = '[#{}:{}]'.format(atom.GetAtomicNum(), atom.GetProp('molAtomMapNumber'))
    else:
        symbol = '[#{}]'.format(atom.GetAtomicNum())

    return symbol

# def convert_atom_to_wildcard(atom):
#     '''Distinguish the difference between aromatic and aliphatic atom. (for supporting information.)'''
#     atom_symbol = atom.GetSymbol()
#     if atom.GetIsAromatic():
#         atom_symbol = atom_symbol.lower()

#     if ':' in atom.GetSmarts():
#         symbol = '[{}:{}]'.format(atom_symbol, atom.GetProp('molAtomMapNumber'))
#     else:
#         symbol = '[{}]'.format(atom_symbol)

#     return symbol


def get_strict_smarts_for_atom(atom):
    '''
    For an RDkit atom object, generate a SMARTS pattern that
    matches the atom as strictly as possible
    
    # 2021/11/16: generalize the atom smarts
    '''
    if ':' in atom.GetSmarts():
        symbol = '[#{}:{}]'.format(atom.GetAtomicNum(), atom.GetProp('molAtomMapNumber'))
    else:
        symbol = '[#{}]'.format(atom.GetAtomicNum())
    
    charge = atom.GetFormalCharge()
    if charge != 0:
        charge_symbol = '+' if (charge >= 0) else '-'
        charge_symbol += '{}'.format(abs(charge))
        if ':' in symbol: 
            symbol = symbol.replace(':', ';{}:'.format(charge_symbol))
        else:
            symbol = symbol.replace(']', ';{}]'.format(charge_symbol))
    
    # if the atom is radical, the number of hydrogen must be included in the template.
    if atom.GetNumRadicalElectrons() != 0:
        num_hydrogen = atom.GetNumExplicitHs()
        if ':' in symbol: 
            symbol = symbol.replace(':', ';H{}:'.format(num_hydrogen))
        else:
            symbol = symbol.replace(']', ';H{}]'.format(num_hydrogen))

    return symbol

# def get_strict_smarts_for_atom(atom):
#     '''
#     Distinguish the difference between aromatic and aliphatic atom. (for supporting information.)
#     # 2024/04/23: Compare.
#     '''
#     atom_symbol = atom.GetSymbol()
#     if atom.GetIsAromatic():
#         atom_symbol = atom_symbol.lower()

#     if ':' in atom.GetSmarts():
#         symbol = '[{}:{}]'.format(atom_symbol, atom.GetProp('molAtomMapNumber'))
#     else:
#         symbol = '[{}]'.format(atom_symbol)
    
#     charge = atom.GetFormalCharge()
#     if charge != 0:
#         charge_symbol = '+' if (charge >= 0) else '-'
#         charge_symbol += '{}'.format(abs(charge))
#         if ':' in symbol: 
#             symbol = symbol.replace(':', ';{}:'.format(charge_symbol))
#         else:
#             symbol = symbol.replace(']', ';{}]'.format(charge_symbol))
    
#     # if the atom is radical, the number of hydrogen must be included in the template.
#     if atom.GetNumRadicalElectrons() != 0:
#         num_hydrogen = atom.GetNumExplicitHs()
#         if ':' in symbol: 
#             symbol = symbol.replace(':', ';H{}:'.format(num_hydrogen))
#         else:
#             symbol = symbol.replace(']', ';H{}]'.format(num_hydrogen))

#     return symbol

def get_fragments_for_changed_atoms(mols, changed_atom_tags, radius=0, 
    category='reactants', expansion=[]):
    '''Given a list of RDKit mols and a list of changed atom tags, this function
    computes the SMILES string of molecular fragments using MolFragmentToSmiles 
    for all changed fragments.

    expansion: atoms added during reactant expansion that should be included and
               generalized in product fragment
    '''
    fragments = ''
    mols_changed = []
    for mol in mols:
        # Initialize list of replacement symbols (updated during expansion)
        symbol_replacements = []

        # Are we looking for special reactive groups? (reactants only)
        if category == 'reactants':
            # groups = get_special_groups(mol)
            groups = []
        else:
            groups = []

        # Build list of atoms to use
        atoms_to_use = []
        for atom in mol.GetAtoms():
            # Check self (only tagged atoms)
            if ':' in atom.GetSmarts():
                if atom.GetSmarts().split(':')[1][:-1] in changed_atom_tags:
                    atoms_to_use.append(atom.GetIdx())
                    symbol = get_strict_smarts_for_atom(atom)
                    if symbol != atom.GetSmarts():
                        symbol_replacements.append((atom.GetIdx(), symbol))
                    continue

        # Fully define leaving groups and this molecule participates?
        if INCLUDE_ALL_UNMAPPED_REACTANT_ATOMS and len(atoms_to_use) > 0:
            if category == 'reactants':
                for atom in mol.GetAtoms():
                    if not atom.HasProp('molAtomMapNumber'):
                        atoms_to_use.append(atom.GetIdx())

        # Check neighbors (any atom)
        for k in range(radius):
            atoms_to_use, symbol_replacements = expand_atoms_to_use(mol, atoms_to_use, 
                groups=groups, symbol_replacements=symbol_replacements)

        if category == 'products':
            # Add extra labels to include (for products only)
            if expansion:
                for atom in mol.GetAtoms():
                    if ':' not in atom.GetSmarts(): continue
                    label = atom.GetSmarts().split(':')[1][:-1]
                    if label in expansion and label not in changed_atom_tags:
                        atoms_to_use.append(atom.GetIdx())
                        # Make the expansion a wildcard
                        symbol_replacements.append((atom.GetIdx(), convert_atom_to_wildcard(atom))) 
                        if VERBOSE: print('expanded label {} to wildcard in products'.format(label))
            
            # Make sure unmapped atoms are included (from products)
            for atom in mol.GetAtoms():
                if not atom.HasProp('molAtomMapNumber'): 
                    atoms_to_use.append(atom.GetIdx())
                    symbol = get_strict_smarts_for_atom(atom)
                    symbol_replacements.append((atom.GetIdx(), symbol))

        # Define new symbols based on symbol_replacements
        symbols = [atom.GetSmarts() for atom in mol.GetAtoms()]
        for (i, symbol) in symbol_replacements:
            symbols[i] = symbol

        if not atoms_to_use: 
            continue
        
        # Keep flipping stereocenters until we are happy...
        # this is a sloppy fix during extraction to achieve consistency
        tetra_consistent = False
        num_tetra_flips = 0
        while not tetra_consistent and num_tetra_flips < 100:
            mol_copy = deepcopy(mol)
            [x.ClearProp('molAtomMapNumber') for x in mol_copy.GetAtoms()]   
            this_fragment = AllChem.MolFragmentToSmiles(mol_copy, atoms_to_use, 
                atomSymbols=symbols, allHsExplicit=True, 
                isomericSmiles=USE_STEREOCHEMISTRY, allBondsExplicit=True)

            # Figure out what atom maps are tetrahedral centers
            # Set isotopes to make sure we're getting the *exact* match we want
            this_fragment_mol = AllChem.MolFromSmarts(this_fragment)
            tetra_map_nums = []
            for atom in this_fragment_mol.GetAtoms():
                if atom.HasProp('molAtomMapNumber'):
                    atom.SetIsotope(int(atom.GetProp('molAtomMapNumber')))
                    if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                        tetra_map_nums.append(atom.GetProp('molAtomMapNumber'))
            map_to_id = {}
            for atom in mol.GetAtoms():
                if atom.HasProp('molAtomMapNumber'):
                    atom.SetIsotope(int(atom.GetProp('molAtomMapNumber')))
                    map_to_id[atom.GetProp('molAtomMapNumber')] = atom.GetIdx()
            
            # Look for matches
            tetra_consistent = True
            all_matched_ids = []
            
            # skip substructure matching if there are a lot of fragments
            # this can help prevent GetSubstructMatches from hanging 
            frag_smi = Chem.MolToSmiles(this_fragment_mol)
            if frag_smi.count('.') > 5:
                break
            
            for matched_ids in mol.GetSubstructMatches(this_fragment_mol, useChirality=True):
                all_matched_ids.extend(matched_ids)
            shuffle(tetra_map_nums)
            for tetra_map_num in tetra_map_nums:
                if VERBOSE: print('Checking consistency of tetrahedral {}'.format(tetra_map_num))
                #print('Using fragment {}'.format(Chem.MolToSmarts(this_fragment_mol, True)))
                if map_to_id[tetra_map_num] not in all_matched_ids:
                    tetra_consistent = False
                    if VERBOSE: print('@@@@@@@@@@@ FRAGMENT DOES NOT MATCH PARENT MOL @@@@@@@@@@@@@@')
                    if VERBOSE: print('@@@@@@@@@@@ FLIPPING CHIRALITY SYMBOL NOW      @@@@@@@@@@@@@@')
                    prevsymbol = symbols[map_to_id[tetra_map_num]]
                    if '@@' in prevsymbol:
                        symbol = prevsymbol.replace('@@', '@')
                    elif '@' in prevsymbol:
                        symbol = prevsymbol.replace('@', '@@')
                    else:
                        raise ValueError('Need to modify symbol of tetra atom without @ or @@??')
                    symbols[map_to_id[tetra_map_num]] = symbol
                    num_tetra_flips += 1
                    # IMPORTANT: only flip one at a time
                    break 

            # Clear isotopes
            for atom in mol.GetAtoms():
                atom.SetIsotope(0)
        
        if not tetra_consistent:
            raise ValueError('Could not find consistent tetrahedral mapping, {} centers'.format(len(tetra_map_nums)))

        fragments += '(' + this_fragment + ').'
        mols_changed.append(Chem.MolToSmiles(clear_mapnum(Chem.MolFromSmiles(Chem.MolToSmiles(mol, True))), True))

    # auxiliary template information: is this an intramolecular reaction or dimerization?
    intra_only = (1 == len(mols_changed))
    dimer_only = (1 == len(set(mols_changed))) and (len(mols_changed) == 2)
    
    return fragments[:-1], intra_only, dimer_only

def extract_from_rxn_smiles(rxn_smiles: str, radius = 0):
    """
    Extract retro template from atom-mapping reaction smiles. 
    This is the generic reaction template extraction.
    """
    reactants = rxn_smiles.split('>')[0]
    products = rxn_smiles.split('>')[-1]
    match_ = Chem.MolFromSmiles(products)
    reactants = mols_from_smiles_list(replace_deuterated(reactants).split('.'))
    products = mols_from_smiles_list(replace_deuterated(products).split('.'))
    
    # if rdkit cant understand molecule, return
    if None in reactants: return None
    if None in products: return None 
    
    # try to sanitize molecules
    try:
        for i in range(len(reactants)):
            reactants[i] = AllChem.RemoveHs(reactants[i]) # *might* not be safe
        for i in range(len(products)):
            products[i] = AllChem.RemoveHs(products[i]) # *might* not be safe
        [Chem.SanitizeMol(mol) for mol in reactants + products] # redundant w/ RemoveHs
        [mol.UpdatePropertyCache() for mol in reactants + products]
    except Exception as e:
        # can't sanitize -> skip
        print(e)
        print('Could not load SMILES or sanitize')
        print(rxn_smiles)
        return None
    
    are_unmapped_product_atoms = False
    extra_reactant_fragment = ''
    for product in products:
        prod_atoms = product.GetAtoms()
        if sum([a.HasProp('molAtomMapNumber') for a in prod_atoms]) < len(prod_atoms):
            if VERBOSE: print('Not all product atoms have atom mapping')
            are_unmapped_product_atoms = True

    if are_unmapped_product_atoms: # add fragment to template
        for product in products:
            prod_atoms = product.GetAtoms()
            # Get unmapped atoms
            unmapped_ids = [
                a.GetIdx() for a in prod_atoms if not a.HasProp('molAtomMapNumber')
            ]
            if len(unmapped_ids) > MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS:
                # Skip this example - too many unmapped product atoms!
                return
            # Define new atom symbols for fragment with atom maps, generalizing fully
            atom_symbols = ['[{}]'.format(a.GetSymbol()) for a in prod_atoms]
            # And bond symbols...
            bond_symbols = ['~' for b in product.GetBonds()]
            if unmapped_ids:
                extra_reactant_fragment += AllChem.MolFragmentToSmiles(
                    product, unmapped_ids, 
                    allHsExplicit = False, isomericSmiles = USE_STEREOCHEMISTRY, 
                    atomSymbols = atom_symbols, bondSymbols = bond_symbols
                ) + '.'
        if extra_reactant_fragment:
            extra_reactant_fragment = extra_reactant_fragment[:-1]
            if VERBOSE: print('extra reactant fragment: {}'.format(extra_reactant_fragment))

        # Consolidate repeated fragments (stoichometry)
        extra_reactant_fragment = '.'.join(sorted(list(set(extra_reactant_fragment.split('.')))))


    if None in reactants + products:
        print('Could not parse all molecules in reaction, skipping')
        print(rxn_smiles)
        return None

    # Calculate changed atoms
    changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
    if err: 
        if VERBOSE:
            print('Could not get changed atoms')
            print(rxn_smiles)
        return
    if not changed_atom_tags:
        if VERBOSE:
            print('No atoms changed?')
            print(rxn_smiles)
        # print('Reaction SMILES: {}'.format(example_doc['RXN_SMILES']))
        return None

    try:
        # Get fragments for reactants
        reactant_fragments, intra_only, dimer_only = get_fragments_for_changed_atoms(reactants, changed_atom_tags, 
            radius = radius, expansion = [], category = 'reactants')
        # Get fragments for products 
        # (WITHOUT matching groups but WITH the addition of reactant fragments)
        product_fragments, _, _  = get_fragments_for_changed_atoms(products, changed_atom_tags, 
            radius = radius, expansion = expand_changed_atom_tags(changed_atom_tags, reactant_fragments),
            category = 'products')
    except ValueError as e:
        if VERBOSE:
            print(e)
            print(rxn_smiles)
        return None

    # Put together and canonicalize (as best as possible)
    rxn_string = '{}>>{}'.format(reactant_fragments, product_fragments)
    rxn_canonical = canonicalize_transform(rxn_string)
    # Change from inter-molecular to intra-molecular 
    rxn_canonical_split = rxn_canonical.split('>>')
    rxn_canonical = rxn_canonical_split[0][1:-1].replace(').(', '.') + \
        '>>' + rxn_canonical_split[1][1:-1].replace(').(', '.')

    reactants_string = rxn_canonical.split('>>')[0]
    products_string  = rxn_canonical.split('>>')[1]

    retro_canonical = products_string + '>>' + reactants_string

    # Load into RDKit
    rxn = AllChem.ReactionFromSmarts(retro_canonical)
    try:
        if rxn.Validate()[1] != 0: 
            print('Could not validate reaction successfully')
            print(rxn_smiles)
            print('retro_canonical: {}'.format(retro_canonical))
            if VERBOSE: input('Pausing...')
            return None
    except Exception as e:
        if VERBOSE:
            print(e)
            print(rxn_smiles)
        return None
    return retro_canonical # retro template

def bond_to_label(bond):
    '''This function takes an RDKit bond and creates a label describing
    the most important attributes
    * This function has been modified. '''
    a1_label = str(bond.GetBeginAtom().GetAtomicNum())
    a2_label = str(bond.GetEndAtom().GetAtomicNum())
    if bond.GetBeginAtom().HasProp('molAtomMapNumber'):
        a1_label += bond.GetBeginAtom().GetProp('molAtomMapNumber')
    if bond.GetEndAtom().HasProp('molAtomMapNumber'):
        a2_label += bond.GetEndAtom().GetProp('molAtomMapNumber')
    atoms = sorted([a1_label, a2_label])

    # modify the second position
    return '{}{}{}'.format(atoms[0], bond.GetBondTypeAsDouble(), atoms[1])

def atoms_are_different(atom1, atom2):
    '''Compares two RDKit atoms based on basic properties'''

    if atom1.GetAtomicNum() != atom2.GetAtomicNum(): return True # must be true for atom mapping
    if atom1.GetTotalNumHs() != atom2.GetTotalNumHs(): return True
    if atom1.GetFormalCharge() != atom2.GetFormalCharge(): return True
    if atom1.GetDegree() != atom2.GetDegree(): return True
    #if atom1.IsInRing() != atom2.IsInRing(): return True # do not want to check this!
    # e.g., in macrocycle formation, don't want the template to include the entire ring structure
    if atom1.GetNumRadicalElectrons() != atom2.GetNumRadicalElectrons(): return True
    if atom1.GetIsAromatic() != atom2.GetIsAromatic(): return True 

    # Check bonds and nearest neighbor identity
    bonds1 = sorted([bond_to_label(bond) for bond in atom1.GetBonds()]) 
    bonds2 = sorted([bond_to_label(bond) for bond in atom2.GetBonds()])
    if bonds1 != bonds2: return True

    return False

if __name__ == "__main__":
    # test_reaction = "[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[NH:8][OH:9].[CH:10]1=[CH:11][CH:12]=[CH:13][CH2:14][CH2:15]1>>[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[O:9][C@H:10]2[CH:11]=[CH:12][C@@H:13]1[CH2:14][CH2:15]2"
    # test_reaction = "CCN([N:3]([CH2:2][CH3:1])[CH2:4][CH:5]=[CH2:6])[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1>>[CH3:1][CH2:2][n:3]1[cH:4][c:5]([CH3:6])[c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]12"
    # template_0 = extract_from_rxn_smiles(test_reaction, radius=0)
    # template_1 = extract_from_rxn_smiles(test_reaction, radius=1)
    # print("test_reaction")
    # print("template radius 0:")
    # print(template_0)
    # print("template radius 1:")
    # print(template_1)


    test_reaction = "[C:3]([c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1)#[CH:10].[CH:11]1=[CH:12][CH:13]2[CH:14]=[CH:15][CH:16]1[CH2:17]2.[O:1]=[C:2]>>[O:1]=[C:2]1[C:3]([c:4]2[cH:5][cH:6][cH:7][cH:8][cH:9]2)=[CH:10][C@H:11]2[C@@H:12]1[C@@H:13]1[CH:14]=[CH:15][C@H:16]2[CH2:17]1"
    template = extract_from_rxn_smiles(test_reaction)
    reac_mols = mols_from_smiles_list(test_reaction.split('>>')[0].split('.'))
    prod_mols = mols_from_smiles_list(test_reaction.split('>>')[1].split('.'))
    for reac_mol in reac_mols:
        find = False
        for reac_atom in reac_mol.GetAtoms():
            if reac_atom.GetAtomMapNum() == 13:
                print('found')
                find = True
                break
        if find:
            break
    
    for prod_mol in prod_mols:
        find = False
        for prod_atom in prod_mol.GetAtoms():
            if prod_atom.GetAtomMapNum() == 13:
                print('found')
                find = True
                break
        if find:
            break
    print(atoms_are_different(reac_atom, prod_atom))
    print(template)