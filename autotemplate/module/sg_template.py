# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 17:51:10 2020

@author: Lung-Yi
"""


# SG_template = [
#     "[C:2]1=[C:4]-[C:5]-[C:9]-[C:8]-[C:3]-1>>[C:2](=[C:3])-[C:4]=[C:5].[C:8]=[C:9]",
#     "[C:2]1=[C:3]-[C:6]2-[C:5]-[C:9]-1-[C:8]=[C:7]-2>>[C:2]#[C:3].[C:5]1-[C:6]=[C:7]-[C:8]=[C:9]-1",
#     "[C:2]1-[C:3]-[C:9]2-[C:8]-[C:12]-1-[C:11]=[C:10]-2>>[C:2](=[C:3]).[C:8]1-[C:9]=[C:10]-[C:11]=[C:12]-1",
#     "[C:2]1=[C:4]-[C:5]-[C:9]=[C:8]-[C:3]-1>>[C:2](=[C:3])-[C:4]=[C:5].[C:8]#[C:9]",
#     "[C:4]1-[C:5]-[C:8]2-[c:9]:[c:10]-[C:11]-1-[c:6]:[c:7]-2>>[C:4]=[C:5].[c:6]1:[c:7]:[c:8]:[c:9]:[c:10]:[c:11]:1",
#     "[C:4]1=[C:5]-[C:8]2-[c:9]:[c:10]-[C:11]-1-[c:6]:[c:7]-2>>[C:4]#[C:5].[c:6]1:[c:7]:[c:8]:[c:9]:[c:10]:[c:11]:1",
#     "[C:2]1=[C:4]-[C:5]-[N:9]-[C:8]-[C:3]-1>>[C:2](=[C:3])-[C:4]=[C:5].[C:8]=[N:9]",
#     "[N:2]1-[C:3]=[C:4]-[C:5]-[C:11]-[C:10]-1>>[N:2]=[C:3]/[C:4]=[C:5].[C:10]=[C:11]",
#     "[C:2]1-[N:3]-[C:9]2-[C:8]-[C:12]-1-[C:11]=[C:10]-2>>[C:2](=[N:3]).[C:8]1-[C:9]=[C:10]-[C:11]=[C:12]-1",
#     "[C:2]1=[C:4]-[N:5]-[C:9]=[C:8]-[C:3]-1>>[C:2](=[C:3])-[C:4]=[N:5].[C:8]#[C:9]",
#     "[N:2]1-[C:3]2-[c:4]:[c:5]-[C:6]-1-[C:12]-[C:13]-2>>[n:2]1:[c:3]:[c:4]:[c:5]:[c:6]:1.[C:12]=[C:13]",
#     "[C:2]1=[C:4]-[C:5]-[#7:9]~[#7:8]-[C:3]-1>>[C:2](=[C:3])-[C:4]=[C:5].[N:8]=[N:9]",
#     "[C:2]1-[S:3]-[C:9]2-[C:8]-[C:12]-1-[C:11]=[C:10]-2>>[C:2](=[S:3]).[C:8]1-[C:9]=[C:10]-[C:11]=[C:12]-1",
#     "[C:2]1=[C:4]-[C:5]-[S:9]-[C:8]-[C:3]-1>>[C:2](=[C:3])-[C:4]=[C:5].[C:8]=[S:9]",
#     "[#6:2]1~[#6:4]-[O:5]-[C:9]-[C:8]-[S:3]-1>>[C:2](=[S:3])-[C:4]=[O:5].[C:8]=[C:9]",
#     "[C:2]1=[C:4]-[C:5]-[C:9]-[C:8]-[O:3]-1>>[C:2](=[O:3])-[C:4]=[C:5].[C:8]=[C:9]",
#     "[C:2]1-[O:9]-[C:10]-[C:8]-[C:6]=[C:4]-1>>[C:2]=[C:4]-[C:6]=[C:8].[O:9]=[C:10]",
#     "[C:2]1=[C:4]-[C:5]-[C:9]=[C:8]-[O:3]-1>>[C:2](=[O:3])-[C:4]=[C:5].[C:8]#[C:9]",
#     "[C:2]1=[C:4]-[C:5]-[O:9]-[N:8](-[C:7]=[O:6])-[C:3]-1>>[C:2](=[C:3])-[C:4]=[C:5].[O:6]=[C:7]-[N:8]-[O:9]",
#     "[C:2]1-[C:7]-[N:8]-[c:9]:[c:10](:[c:11])-[C:4]-1-[O:5]>>[C:2]-[C:4]=[O:5].[C:7]=[N:8]-[c:9]:[c:10]:[c:11]",
#     "[N:3]1-[O:4]-[C:8]-[C:7]=[C:6]-[C:10]-1>>[N:3]=[O:4].[C:10]=[C:6]-[C:7]=[C:8]",
#     ]

SG_template_Diels_Alder = [
    "[*:3]1-[*:4]~[*:5]-[*:6]-[*:1]~[*:2]-1>>[*:3]=[*:4]-[*:5]=[*:6].[*:1]=[*:2]",
    "[*:3]1-[*:4]~[*:5]-[*:6]-[*:1]~[*:2]-1>>[*:3]=[*:4]-[*:5]=[*:6].[*:1]#[*:2]",
    "[C:2]1=[C:4]-[C:5]-[O:9]-[N:8](-[C:7]=[O:6])-[C:3]-1>>[C:2](=[C:3])-[C:4]=[C:5].[O:6]=[C:7]-[N:8]-[O:9]",
    "[C:2]1-[C:4](-[O:5])-[c:7](:[c:6]):[c:8]-[N:9]-[C:10]-1>>[C:2]-[C:4]=[O:5].[c:6]:[c:7]:[c:8]-[N:9]=[C:10]",
    "[N:2]1-[C:3]=[C:4]-[C:5]-[C:11](-[C:10])-[C:12]-1=[O:13]>>[N:2]=[C:3]-[C:4]=[C:5].[C:10]=[C:11]-[C:12]=[O:13]",
    "[C:12]1-[C:2]2-[C:3]=[C:4]-[C:5](-[C:6]-[C:1]-2)-[C:13]-1>>[C:1]=[C:2]-[C:3]=[C:4]-[C:5]=[C:6].[C:12]=[C:13]"
    ]

SG_template_Reductive_Amination = [
    "[*:3]-[#7:4]-[#6:1]>>O=[#6:1].[*:3]-[#7:4]",
    "[*:7]-[N:8](-[C:2])-[C:5]>>[C:2]=O.[C:5]=O.[*:7]-[N:8]",
    "[#6:2]-[#7:1]-[#6:8]-[#7:4]-[#6:5]>>[#7:1]-[#6:2].[#7:4]-[#6:5].O=[#6:8]"
    ]

SG_template_Hydrolysis = [
    "[#6:4]-[*:1]-[*:2]=[O:3]>>[O:3]=[*:2]-O.[#6:4]-[*:1]",
    "[C:2]-[O,N:1]-[C:6](-[C:5])-[O,N:4]-[C:3]>>[O,N:1]-[C:2].[C:3]-[O,N:4].[C:5]-[C:6](=O)",
    "[C:2]-[O,N:1]-[C:6](-[C:5])-[O,N:4]-[C:3]>>[O,N:1]-[C:2]-[C:3]-[O,N:4].[C:5]-[C:6](=O)",
    #"[C:2]-[O:1]-[C:6](-[C:5])-[O:4]-[C:3]>>[O:1]-[C:2]~[C:3]-[O:4].[C:5]-[C:6](=O)",
    "[N:4]=[C:1](-[*:2])-[*:3]>>[*:2]-[C:1](=O)-[*:3].[N:4]",
    "[N:1]#[C:2]>>[N:1]-[C:2]=O",
    "[N:3]#[C:1]>>[C:1](=O)-O.[N:3]",
    "[Si:3]-[O:1]-[C:2]>>[C:2]-[O:1].[Si:3]-O",
    "[C:3](=[O:1])-[N:2]>>[N:2].[C:3](=[O:1])-O",
    "[#6:1]-[F,Cl,Br,I:2]>>[#6:1]-O.[F,Cl,Br,I:2]",
    "[O:1]=[C:2]-[C:4]=[C:5]>>[O:1]-[C:2]-[C:4]-[C:5]-O",
    "[C:2](=[N:3])-[N:5]>>O=[C:2]-[N:3].[N:5]"
    #"[C:2]-[O:1]-[C:6](-[C:5])-[O:4]-[C:3]>>[O:1]-[C:2].[C:3]-[O:4].[C:5]-[C:6](=O)"
    #"[C:2]-[O:1]-[C:6](-[C:7])(-[C:5])-[O:4]-[C:3]>>[O:1]-[C:2].[C:3]-[O:4].[C:5]-[C:6](=O)-[C:7]"
    ]

SG_template_Hydrogenation = [
    "[C:1]=[*:2]>>[C:1]-[*:2]",
    "[C:1]#[*:2]>>[C:1]-[*:2]",
    "[O:3]=[N+;H0;D3:1](-[O-:4])-[#6:2]>>[NH2;D1;+0:1]-[#6:2].[O:3].[O;+0:4]",
    "[C:1]#[*:2]>>[C:1]=[*:2]",
    "[c,n,o,s:1]1:[c,n,o,s:2]:[c,n,o,s:3]:[c,n,o,s:4]:[c,n,o,s:5]:[c,n,o,s:6]:1>>[C,N,O,S:1]1-[C,N,O,S:2]-[C,N,O,S:3]-[C,N,O,S:4]-[C,N,O,S:5]-[C,N,O,S:6]-1",
    "[*:1]-[Cl,Br,F:2]>>[*:1].[Cl,Br,F:2]",
    "[C:1](=[O:2])-[O:3]>>[C:1]-[O:3].[O:2]",
    "[c,n,o,s:1]1:[c,n,o,s:2]:[c,n,o,s:3]:[c,n,o,s:4]:[c,n,o,s:5]:1>>[C,N,O,S:1]1-[C,N,O,S:2]-[C,N,O,S:3]-[C,N,O,S:4]-[C,N,O,S:5]-1",
    "[#6:1]-[#6:2](=[O:4])-[#6:3]>>[#6:1]-[#6:2]-[#6:3]",
    "[O:3]-[N:1]=[C:2]>>[C:2]-[N:1]",
    "[O:1]=[C:2]-[O:3]-[C:4]>>[O:1]=[C:2]-[O:3].[C:4]",
    "[C:2]1-[O:4]-[C:5]-1>>[C:2]-[C:5]-[O:4]",
    "[C:2]-[O:4]-[C:5]>>[C:2].[C:5]-[O:4]",
    "[N;H0;D2;+0:1]=[O:2]>>[NH2;D1;+0:1]",
    "[C:2]1=[C:3]-[O:5]-[C:6](=[O:7])-[C:8]-1>>[C:2](-[C:8]-[C:6](=[O:7])-[O:5])-[C:3]",
    "[C:2]1-[C:3]-[C:4]-1>>[C:2]-[C:3]-[C:4]",
    "[c,n,o,s:1]1:[c,n,o,s:2]:[c,n,o,s:3]:[c,n,o,s:4]:[c,n,o,s:5]:[c,n,o,s:6]:1>>[C,N,O,S:1]1-[C,N,O,S:2]-[C,N,O,S:3]-[C,N,O,S:4]-[C,N,O,S:5]:[C,N,O,S:6]-1",
    "[c,n,o,s:1]1:[c,n,o,s:2]:[c,n,o,s:3]:[c,n,o,s:4]:[c,n,o,s:5]:1>>[C,N,O,S:1]1-[C,N,O,S:2]-[C,N,O,S:3]-[C,N,O,S:4]:[C,N,O,S:5]-1",
    "[O:1]-[O:2]>>[O:1].[O:2]"
    ]

"""
start testing
"""
# def rdchiralRunText_modified(reaction_smarts, reactant_smiles):
#     rxn = rdchiralReaction(reaction_smarts)
#     rxn.reset()
#     reactants = rdchiralReactants(reactant_smiles)
#     # Run naive RDKit on ACHIRAL version of molecules
#     outcomes = rxn.rxn.RunReactants((reactants.reactants_achiral,))
#     # mol = Chem.MolFromSmiles(reactant_smiles)
#     # Chem.rdmolops.RemoveStereochemistry(mol)
#     # outcomes = rxn.rxn.RunReactants((mol,)) # for tesret
#     smiles_list = []
#     for outcome in outcomes:
#         ###############################################################################
#         # Look for new atoms in products that were not in 
#         # reactants (e.g., LGs for a retro reaction)
#         changed_index = []
        
#         unmapped = sum([z.GetNumAtoms() for z in outcome])
#         for m in outcome:
#             try:
#                 for a in m.GetAtoms():
#                     # Assign map number to outcome based on react_atom_idx
#                     if a.HasProp('react_atom_idx'):
#                         num = reactants.idx_to_mapnum(int(a.GetProp('react_atom_idx')))
#                         if a.GetAtomMapNum() != num: 
#                             changed_index.append(num)
#                             a.UpdatePropertyCache(strict=False)
#                             try:
#                                 if a.GetTotalValence() != atom_valence_dict[a.GetSymbol()]:
#                                     a.SetNumExplicitHs(CalculateNumHs(a))
#                                     a.UpdatePropertyCache(strict=True)
#                             except:
#                                 pass
#                             # UpdateAtomInfo(a)
#                         a.SetAtomMapNum(num)
#                     if not a.GetAtomMapNum():
#                         a.SetAtomMapNum(unmapped)
#                         unmapped -= 1
#             except:
#                 continue
#         smiles = ""
#         for mol in outcome:
#             smiles += Chem.MolToSmiles(mol) + "."
        
#         smiles_list.append(smiles.rstrip("."))
#     return smiles_list