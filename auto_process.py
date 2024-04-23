import argparse
import os
from tqdm import tqdm
from rdkit import Chem
import matplotlib.pyplot as plt
import pickle

from collections import Counter
from autotemplate.extract_utils import extract_from_rxn_smiles, canon_remap
from autotemplate.run_utils import ReassignMapping, rdchiralRunText_modified, RemoveReagent, clearIsotope, build_removed_reagents_dict
from autotemplate.graph_utils import mapping_for_gold_multiple_smiles, find_unique_templates_dict

def count_change(rxn_smiles):
    rea, prod = rxn_smiles.split('>>')
    rea_count = rea.count('.') + 1
    prod_count = prod.count('.') + 1
    return '{}>>{}'.format(rea_count, prod_count)

def IsInvalidSmiles(smiles):
    """ Check whether there exists repeated atom mapping number in the SMILES.
        Or the rdkit of molecule is None"""
    mol = Chem.MolFromSmiles(smiles)
    if mol == None:
        return True
    atom_map = []
    for atom in mol.GetAtoms():
        mapnum = atom.GetProp('molAtomMapNumber')
        if mapnum in atom_map:
            print(atom_map + [mapnum])
            return True
        atom_map.append(mapnum)
    else:
        return False

def remap_all_templates(reactants, product, all_templates, retro = True, return_ChiralStereo=True):
    """
    This function is retro template, because we can quickly examine if there is secondary amine functional group
    in product. Some reactions are devoid of reactant, so we need to append it on the reactant site.
    """
    r_gold_smiles = reactants
    r_check, r_numatom = canon_remap(reactants, return_NumAtom=True)
    product, p_numatom = ReassignMapping(product, return_NumAtom=True, iso=return_ChiralStereo)

    # if t == None:
    for template in all_templates:
        outcomes = remap_one_template(template, product)
        if not outcomes: continue

        for outcome in outcomes:
            if IsInvalidSmiles(outcome): continue
            r = canon_remap(outcome)
            if r == None: continue
            r_copy = r.split('.').copy()
            if not ([i for i in r_check.split('.') if not i in r_copy or r_copy.remove(i)]):
            #if not (False in [s in r.split('.') for s in r_check.split('.')]):
            #all Reaxys reactant records must lie in template-outcome
                if return_ChiralStereo:
                    add = r.split('.')
                    for s in r_check.split('.'): add.remove(s)
                    if add: r_gold_smiles = reactants + '.' + '.'.join(add)
                    gold_outcome = mapping_for_gold_multiple_smiles(r_gold_smiles, outcome)
                    return gold_outcome+">>"+product, template
                else:
                    return outcome+">>"+product, template
    return None, None

def remap_one_template(template, target_smiles):
    try:
        outcomes = rdchiralRunText_modified(template, target_smiles)
        return outcomes
    except:
        return

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

def check_reconstruct(template, rxn_smiles, retro = True):
    """Check whether the template obtained from the function "extract_from_rxn_smiles" is True for the original rxn_smiles. """
    reactants, products = rxn_smiles.split('>>')
    gold_reac = set(canon_remap(reactants).split('.'))
    
    if retro:
        reac_list = remap_one_template(template, products)
        if not reac_list: return False
        for reac in reac_list:
            reac = canon_remap(reac)
            if reac:
                check_reac = set(reac.split('.'))
            else:
                continue
            if gold_reac.issubset(check_reac):
                return True
        return False            
    else:
        pass
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file',type=str,
                        default='./data/PausonKhandReaction/MappingResult_PausonKhandReaction.txt')
    parser.add_argument('--minimum_percent',type=float,
                        default=0.01)
    parser.add_argument('--radius',type=int,
                        default=0  )
    parser.add_argument('--threshold',type=int,
                        default=5)
    parser.add_argument('--skip_topn_templates',type=int,
                        default=3)
    parser.add_argument('--retain_reagents',type=bool,
                        default=True)
    # parser.add_argument('--verbose', default=False, type=lambda x: (str(x).lower() in ['true','1', 'yes'])) # Whether to print failed prediction information
    args = parser.parse_args()

    f = open(args.input_file, 'r')
    input_file = f.readlines()
    f.close()
    
    """ Remove reagent. 
    Remove 
    (1) same molecule appear in both reactant and product sites. 
    (2) unmapped molecules in both reactant and product sites. """
    print('Removing reagent...')
    input_file = [[clearIsotope(line.split('\t')[0]), line.split('\t')[1].strip('\n')] for line in input_file]
    removed_reagents_list = []
    for i, (rxn_smiles, reaction_id) in enumerate(input_file):
        if args.retain_reagents:
            processed_rxn_smiles, removed_reagents = RemoveReagent(rxn_smiles, retain_reagents=args.retain_reagents)
            removed_reagents_list.append(removed_reagents)
        else:
            processed_rxn_smiles = RemoveReagent(rxn_smiles, retain_reagents=args.retain_reagents)
        input_file[i][0] = processed_rxn_smiles

    if args.retain_reagents:
        id2reagent_dict = build_removed_reagents_dict(list(zip(*input_file))[1], removed_reagents_list)
        print(id2reagent_dict)
    
    """ Analyze the reaction change """
            
    count_dict = Counter([count_change(rxn_smiles) for rxn_smiles, _ in input_file])
    joint = list(count_dict.items())
    joint = sorted(joint, key = lambda x: x[0])
    x, y =zip(*joint)
    plt.figure(dpi=500)
    plt.bar(x,y)
    plt.xticks(rotation=45)
    plt.savefig(os.path.join(os.path.dirname(args.input_file), 'reaction.png'))
    
    """ Take redictive amination as example reaction.
        We find that it is suitable to use retro template to preprocess this reaciton..."""
    """ The extracted templates may share the exactly same reaction change but with different reaction SMARTS expression.
        Therefore, we need to standarize these templates. """
    
    USE_RETRO_TEMPLATE = True
    
    print('Extracting templates...')
    templates = []
    data = []
    for i in tqdm(range(len(input_file))):
        rxn_smiles, reaxys_id = input_file[i]
        template = extract_from_rxn_smiles(rxn_smiles, radius = args.radius)

        # TODO: avoid using the wrong atom-mapped reaction template
        if check_reconstruct(template, rxn_smiles, retro = True):
            templates.append(template)
        #     data.append({'rxn_smiles': rxn_smiles, 'template': template, 'reaxys id': reaxys_id})
        # else:
        data.append({'rxn_smiles': rxn_smiles, 'template': None, 'reaxys id': reaxys_id})
        
    templates = Counter(templates)
    templates, changed_records = find_unique_templates_dict(templates)
    # process the data of using what template:
    for i in range(len(data)):
        old_used_template = data[i]['template']
        if old_used_template:
            data[i]['template'] = changed_records[old_used_template]
    
    final_templates = list(templates.items())
    final_templates = sorted(final_templates, key = lambda x: x[1], reverse = True)
    _, num = zip(*final_templates)
    num = list(num)
    
    if args.threshold:
        threshold = args.threshold
    else:
        threshold = int(len(data)*args.minimum_percent)
    final_templates = [ x for x in final_templates if (x[1] > threshold) and x[0] != None]
    all_templates = list(zip(*final_templates))[0]

    
    """ Start to apply the final templates on the reaction with uncommon template."""
    print('Applying final templates on the wrong reaction...')
    changed_data = []
    for i in tqdm(range(len(data))):
        rxn_smiles = data[i]['rxn_smiles']
        its_template = data[i]['template']
        if its_template in all_templates[:args.skip_topn_templates]:
            new_rxn_smiles = give_mapnum_for_unmapped_reactant_atoms(rxn_smiles)
            data[i].update({'final rxn_smiles': new_rxn_smiles})
            continue # we do not modify the reaction with the correct tmeplate.
        else:
            reactants, product = rxn_smiles.split('>>')
            if '.' in product:
                for subproduct in product.split('.'): # apply templates on single product reaction
                    templated_rxn, used_template = remap_all_templates(reactants, subproduct, all_templates, retro = USE_RETRO_TEMPLATE, return_ChiralStereo = True)
                    if templated_rxn:
                        data[i].update({'final rxn_smiles': templated_rxn, 'template': used_template})
                        changed_data.append(data[i])
                        break
                else:
                    data[i].update({'final rxn_smiles': templated_rxn, 'template': used_template})
                    changed_data.append(data[i])
            else:
                templated_rxn, used_template = remap_all_templates(reactants, product, all_templates, retro = USE_RETRO_TEMPLATE, return_ChiralStereo = True)
                data[i].update({'final rxn_smiles': templated_rxn, 'template': used_template})
                changed_data.append(data[i])
    
    """ Write template processed data. """
    g = open(args.input_file + '.processed', 'w')
    h = open(args.input_file + '.failed', 'w')
    for i in range(len(data)):
        rxn_smiles = data[i]['final rxn_smiles']
        if rxn_smiles:
            if args.retain_reagents:
                reagents = id2reagent_dict[data[i]['reaxys id']]
                if reagents:
                    rxn_smiles += "."+reagents
            g.write(rxn_smiles + '\t' + data[i]['reaxys id'] + '\n')
        else:
            no_processed_smiles = data[i]['rxn_smiles']
            if args.retain_reagents:
                reagents = id2reagent_dict[data[i]['reaxys id']]
                if reagents:
                    no_processed_smiles += "."+reagents
            h.write(no_processed_smiles + '\t' + data[i]['reaxys id'] + '\n')
    g.close()
    h.close()
    
    with open(os.path.join(os.path.dirname(args.input_file), 'processed_data.pkl'), 'wb') as f:
        pickle.dump(data, f)
    
    """ Write used templates. """
    template_output = os.path.join(os.path.dirname(args.input_file), 'all_templates_used.csv')
    z = open(template_output, 'w')
    z.write("template"+','+"count"+'\n')
    for template, count in final_templates:
        z.write(template+','+str(count)+'\n')
    z.close()
    
