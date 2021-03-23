import os
import copy

PATH_TO_BASE_FOLDER = '/home/fol007/PhD_Project/Template_Based_Docking_Project/cross_free_docking_cluster_run'
PATH_TO_REFERENCE_LIGANDS_FOLDER = f'{PATH_TO_BASE_FOLDER}/reference_ligands'
PATH_TO_DONE_FOLDER = f'{PATH_TO_BASE_FOLDER}/aligned/done'


def create_initial_reference_dictionary(path_to_reference_ligands_folder=PATH_TO_REFERENCE_LIGANDS_FOLDER):
    '''organize dictionary based on the reference ligands'''
    reference_dictionary = {}
    uniprot_ids = os.listdir(path_to_reference_ligands_folder)
    for uniprot_id in uniprot_ids:
        reference_dictionary[uniprot_id] = []
        for references in os.listdir(path_to_reference_ligands_folder + '/' + uniprot_id):
            reference_dictionary[uniprot_id] += [references.split('.')[0]]
    return reference_dictionary


inputs = create_initial_reference_dictionary()
inputs_copy = copy.deepcopy(inputs)
output = create_initial_reference_dictionary(path_to_reference_ligands_folder=PATH_TO_DONE_FOLDER)

for uniprot_id in output.keys():
    # if they are the same length they everything is done regarding that target
    if len(inputs[uniprot_id]) == len(output[uniprot_id]):
        inputs_copy.pop(uniprot_id)
        continue
    # if some are missing, remove the ones that were already done
    for already_done in output[uniprot_id]:
        inputs_copy[uniprot_id].remove(already_done)

if len(inputs_copy) == 0:
    print('All tasks completed!')
else:
    for target in inputs_copy:
        total = len(inputs[target])
        left_to_do = len(inputs_copy[target])
        print(' '.join([target, 'has', str(left_to_do), 'out of', str(total), 'left to do!']))
