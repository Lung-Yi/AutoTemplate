import argparse
import os
from tqdm import tqdm
from rdkit import Chem
import matplotlib.pyplot as plt
import pickle
import multiprocessing

from collections import Counter
from autotemplate.extract_utils import extract_from_rxn_smiles, canon_remap
from autotemplate.run_utils import ReassignMapping, rdchiralRunText_modified, RemoveReagent, clearIsotope, build_removed_reagents_dict
from autotemplate.graph_utils import mapping_for_gold_multiple_smiles, find_unique_templates_dict, give_mapnum_for_unmapped_reactant_atoms, balance_rxn_smiles

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
    parser.add_argument('--ncpus',type=int,
                        default=1)
    parser.add_argument('--retain_reagents', action='store_true')
    parser.add_argument('--save_threshold', action='store_true')
    parser.add_argument('--balance_product', action='store_true')
    args = parser.parse_args()

    f = open(args.input_file, 'r')
    input_file = f.readlines()
    f.close()
    
    """ Remove reagent. 
    Remove 
    (1) same molecule appear in both reactant and product sites. 
    (2) unmapped molecules in both reactant and product sites. """
    print('Removing reagent...')
    def subfunction_clearIsotope(line):
        return [clearIsotope(line.split('\t')[0]), line.split('\t')[1].strip('\n')]
    with multiprocessing.Pool(processes=args.ncpus) as pool:
        input_file = pool.map(subfunction_clearIsotope, input_file)

    if args.retain_reagents:
        def subfunction_remove_reagents(rxn_smiles):
            processed_rxn_smiles, removed_reagents = RemoveReagent(rxn_smiles, retain_reagents=True)
            return processed_rxn_smiles, removed_reagents
    else:
        def subfunction_remove_reagents(rxn_smiles):
            processed_rxn_smiles = RemoveReagent(rxn_smiles, retain_reagents=False)
            return processed_rxn_smiles
        
    all_rxn_smiles = list(zip(*input_file))[0]
    with multiprocessing.Pool(processes=args.ncpus) as pool:
        processed_results = pool.map(subfunction_remove_reagents, all_rxn_smiles)
    
    # update input file:
    for i in range(len(input_file)):
        if args.retain_reagents:
            input_file[i][0] = processed_results[i][0]
        else:
            input_file[i][0] = processed_results[i]

    if args.retain_reagents:
        id2reagent_dict = build_removed_reagents_dict(list(zip(*input_file))[1], list(zip(*processed_results))[1])
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
    
    """ Take reductive amination as example reaction.
        We find that it is suitable to use retro template to preprocess this reaciton..."""
    """ The extracted templates may share the exactly same reaction change but with different reaction SMARTS expression.
        Therefore, we need to standarize these templates. """
    
    USE_RETRO_TEMPLATE = True
    
    print('Extracting templates...')
    def subfunction(rxn_smiles):
        """This is a child function for multiprocessing."""
        template = extract_from_rxn_smiles(rxn_smiles, radius = args.radius)
        if check_reconstruct(template, rxn_smiles, retro = True):
            return template
        else:
            return None
    
    rxn_smiles_list = []
    templates = []
    for i in range(len(input_file)):
        rxn_smiles, reaxys_id = input_file[i]
        rxn_smiles_list.append(rxn_smiles)

    with multiprocessing.Pool(processes=args.ncpus) as pool:
        templates = pool.map(subfunction, rxn_smiles_list)
    templates = [i for i in templates if i is not None]

    data = []
    for i in range(len(input_file)):
        rxn_smiles, reaxys_id = input_file[i]
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
    print(f"Number of generic templates: {len(all_templates)}")

    
    """ Start to apply the final templates on the reaction with uncommon template."""
    print('Applying final templates on the wrong reaction...')
    def sub_function_applying_generic_templates(rxn_smiles):
        """ This the child function for multiprocessing. It is used to iteratively applying generic templates on the product."""
        reactants, product = rxn_smiles.split('>>')
        if '.' in product:
            for subproduct in product.split('.'): # apply templates on single product reaction
                templated_rxn, used_template = remap_all_templates(reactants, subproduct, all_templates, retro = USE_RETRO_TEMPLATE, return_ChiralStereo = True)
                if templated_rxn:
                    return templated_rxn, used_template
            else:
                return templated_rxn, used_template
        else:
            templated_rxn, used_template = remap_all_templates(reactants, product, all_templates, retro = USE_RETRO_TEMPLATE, return_ChiralStereo = True)
            return templated_rxn, used_template
        
    
    update_index = [] # the list index of the data that needs to be updated by applying generic tepmlates.
    update_rxn_smiles_list = []
    for i in range(len(data)):
        rxn_smiles = data[i]['rxn_smiles']
        its_template = data[i]['template']
        if its_template in all_templates[:args.skip_topn_templates]:
            new_rxn_smiles = give_mapnum_for_unmapped_reactant_atoms(rxn_smiles)
            data[i].update({'final rxn_smiles': new_rxn_smiles})
            continue # we do not modify the reaction with the correct tmeplate. (Assume the top-n templates are correct because of their high frequency.)
        else:
            update_index.append(i)
            update_rxn_smiles_list.append(rxn_smiles)
    
    with multiprocessing.Pool(processes=args.ncpus) as pool:
        results_list = pool.map(sub_function_applying_generic_templates, update_rxn_smiles_list)
    
    for index, result in zip(update_index, results_list):
        templated_rxn, used_template = result
        data[index].update({'final rxn_smiles': templated_rxn, 'template': used_template})
    
    """ Write template processed data. """
    if args.save_threshold:
        save_dir = os.path.join(os.path.dirname(args.input_file), f"threshold_{args.threshold}") + "/"
        os.makedirs(save_dir, exist_ok=True)
        head, tail = os.path.split(args.input_file)
        g = open(save_dir + tail + '.processed', 'w')
        h = open(save_dir + tail + '.failed', 'w')
    else:
        g = open(args.input_file + '.processed', 'w')
        h = open(args.input_file + '.failed', 'w')
    
    for i in range(len(data)):
        rxn_smiles = data[i]['final rxn_smiles']
        if rxn_smiles:
            if args.balance_product:
                rxn_smiles = balance_rxn_smiles(rxn_smiles)
            if args.retain_reagents:
                reagents = id2reagent_dict[data[i]['reaxys id']]
                if reagents:
                    rxn_smiles = reagents + "." + rxn_smiles
            g.write(rxn_smiles + '\t' + data[i]['reaxys id'] + '\n')
        else:
            no_processed_smiles = data[i]['rxn_smiles']
            if args.retain_reagents:
                reagents = id2reagent_dict[data[i]['reaxys id']]
                if reagents:
                    rxn_smiles = reagents + "." + rxn_smiles
            h.write(no_processed_smiles + '\t' + data[i]['reaxys id'] + '\n')
    g.close()
    h.close()
    
    if args.save_threshold:
        with open(save_dir + 'processed_data.pkl', 'wb') as f:
            pickle.dump(data, f)
    else:
        with open(os.path.join(os.path.dirname(args.input_file), 'processed_data.pkl'), 'wb') as f:
            pickle.dump(data, f)
    
    """ Write used templates. """
    if args.save_threshold:
        template_output = save_dir + 'all_templates_used.csv'
    else:
        template_output = os.path.join(os.path.dirname(args.input_file), 'all_templates_used.csv')
    z = open(template_output, 'w')
    z.write("template"+','+"count"+'\n')
    for template, count in final_templates:
        z.write(template+','+str(count)+'\n')
    z.close()
    
