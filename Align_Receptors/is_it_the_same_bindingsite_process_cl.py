
# import necessary modules
import os
import multiprocessing as mp
import itertools
import numpy as np
import pandas as pd
from pymol import cmd

# define an output queue
output = mp.Queue()

PATH_TO_BASE_FOLDER = '/home/fol007/PycharmProjects/ChEMBL_plus_BindingMOAD/BindingMOAD_AstexDiverseSet_Simplified'#'/BindingMOAD_AstexDiverseSet_Simplified'
PATH_TO_PDB_FOLDER = f'{PATH_TO_BASE_FOLDER}/pdb_files'
PATH_TO_REFERENCE_LIGANDS_FOLDER = f'{PATH_TO_BASE_FOLDER}/reference_ligands'
PATH_TO_GUIDE = f'{PATH_TO_BASE_FOLDER}/../../Analysis_of_Docking/data/rmsd_values_base.csv'

def main():
    # get all uniprot ids
    uniprots_in_guide = pd.read_csv(PATH_TO_GUIDE)['uniprot_id'].drop_duplicates().values

    reference_dictionary = create_reference_dictionary(list_of_uniprots=uniprots_in_guide)

    dictionary_df = {'uniprot_id':[], 'ref1':[], 'ref2':[]}

    for uniprot_id in reference_dictionary:
        for ref1, ref2 in reference_dictionary[uniprot_id]:
            dictionary_df['uniprot_id'] += [uniprot_id]
            dictionary_df['ref1'] += [ref1.split('.')[0]]
            dictionary_df['ref2'] += [ref2.split('.')[0]]

    dictionary_df = pd.DataFrame(dictionary_df)

    dictionary_df['rmsd'] = None
    dictionary_df['distance'] = None
    dictionary_df['error'] = None

    processes = []
    i = 0
    for index, uniprot_id, ref1, ref2 in dictionary_df[['uniprot_id', 'ref1', 'ref2']].itertuples():

        path_to_ligang_ref1 = PATH_TO_REFERENCE_LIGANDS_FOLDER + '/' + uniprot_id + '/' + ref1 + '.sdf'
        path_to_ligand_ref2 = PATH_TO_REFERENCE_LIGANDS_FOLDER + '/' + uniprot_id + '/' + ref2 + '.sdf'
        path_to_pdb_ref1 = PATH_TO_PDB_FOLDER + '/' + uniprot_id + '/' + ref1.split('_')[0] + '.pdb'
        path_to_pdb_ref2 = PATH_TO_PDB_FOLDER + '/' + uniprot_id + '/' + ref2.split('_')[0] + '.pdb'

        processes += [mp.Process(target=run_in_parallel, args=(index, path_to_ligang_ref1, path_to_ligand_ref2, path_to_pdb_ref1, path_to_pdb_ref2, output))]

        i += 1
        if i ==50:
            break

    # Run processes
    for p in processes:
        p.start()

    # Exit the completed processes
    for p in processes:
        p.join()

    # Get process results from the output queue
    results = [output.get() for p in processes]

    for index, rmsd, distance, error in results:
        dictionary_df.at[index, 'rmsd'] = rmsd
        dictionary_df.at[index, 'distance'] = distance
        dictionary_df.at[index, 'error'] = error

    dictionary_df.to_csv('dictionary.csv', index=False)

def run_in_parallel(index, path_to_ligang_ref1, path_to_ligand_ref2, path_to_pdb_ref1, path_to_pdb_ref2, output):
    try:
        rmsd, distance = align(path_to_ligang_ref1, path_to_ligand_ref2, path_to_pdb_ref1, path_to_pdb_ref2)
        error = None
    except Exception as e:
        rmsd = None
        distance = None
        error = e
    output.put((index, rmsd, distance, error))

def align(path_to_ligang_ref1, path_to_ligang_ref2, path_to_pdb_ref1, path_to_pdb_ref2):

    cmd.load(path_to_pdb_ref1, object='2gtn')
    cmd.load(path_to_pdb_ref2, object='3p7b')
    cmd.load(path_to_ligang_ref1, object='lie')
    cmd.load(path_to_ligang_ref2, object='p7b')

    cmd.alter(selection='lie', expression='resn="lie"')
    cmd.alter(selection='p7b', expression='resn="p7b"')

    cmd.extract('2gtn_lie', selection='(2gtn, lie)')
    cmd.extract('3p7b_p7b', selection='(3p7b, p7b)')

    cmd.select(name='lie_receptor', selection='br. 2gtn_lie within 10 of /2gtn_lie///lie')
    cmd.select(name='p7b_receptor', selection='br. 3p7b_p7b within 10 of /3p7b_p7b///p7b')

    params = cmd.align('lie_receptor', 'p7b_receptor')
    rmsd = params[0]

    com1 = cmd.centerofmass(selection='/3p7b_p7b///p7b')

    com2 = cmd.centerofmass(selection='/2gtn_lie///lie')

    d = np.array(com1) - np.array(com2)
    com1_to_com2 = np.sqrt(np.dot(d, d))

    cmd.reinitialize()

    return rmsd, com1_to_com2

def create_reference_dictionary(list_of_uniprots, path_to_reference_ligands_folder=PATH_TO_REFERENCE_LIGANDS_FOLDER):
    '''organize dictionary based on the reference ligands'''
    reference_dictionary = {}
    uniprot_ids = os.listdir(path_to_reference_ligands_folder)
    for uniprot_id in uniprot_ids:
        if uniprot_id in list_of_uniprots:
            reference_dictionary[uniprot_id] = []
            for references in os.listdir(path_to_reference_ligands_folder + '/' + uniprot_id):
                reference_dictionary[uniprot_id] += [references]
            reference_dictionary[uniprot_id] = list(itertools.combinations(reference_dictionary[uniprot_id], 2))
    return reference_dictionary

if __name__ == "__main__":

    main()
