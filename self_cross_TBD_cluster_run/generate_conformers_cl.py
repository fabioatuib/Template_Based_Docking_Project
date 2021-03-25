
# import necessary modules
import time
import os
import sys
import numpy as np
import multiprocessing
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem

# define constants
MAX_CONFS = 1000

PATH_TO_BASE_FOLDER = '/run/user/121976/gvfs/sftp:host=biomed3100111/Home/siv32/fol007/new_cross_docking_protocol/Self_and_Cross_TBD'
PATH_TO_REFERENCE_LIGANDS_FOLDER = f'{PATH_TO_BASE_FOLDER}/reference_ligands'
PATH_TO_SMILES_TO_GENERATE_CONFORMER = f'{PATH_TO_BASE_FOLDER}/smiles_to_genconformer'
PATH_TO_CONFORMERS_TO_DOCK = f'{PATH_TO_BASE_FOLDER}/conformers_to_dock'
PATH_TO_LOG = f'{PATH_TO_BASE_FOLDER}/logs_conformer_generation'

def main():

    initial_reference_dictionary = create_initial_reference_dictionary()
    final_reference_dictionary = create_final_reference_dictionary()

    already_done = []
    for uniprot_id in final_reference_dictionary.keys():
        if len(initial_reference_dictionary[uniprot_id]) == len(final_reference_dictionary[uniprot_id]):
            already_done += [uniprot_id]

    print(initial_reference_dictionary)
    print(already_done)

    uniprot_count = [(key, len(initial_reference_dictionary[key])) for key in initial_reference_dictionary.keys()]
    uniprot_count.sort(key=lambda x: x[1])

    for uniprot_id, count in uniprot_count:
        if uniprot_id in already_done:
            continue
        sub_main(uniprot_id, count)

def sub_main(uniprot_id, number_of_ligands):
    t = time.time()

    log = open(PATH_TO_LOG + '/' + uniprot_id + '.txt', 'w')
    log.write('uniprot_id: ' + uniprot_id + '\n')
    log.write('number_of_ligands: ' + str(number_of_ligands) + '\n' + '' + '\n')
    # 1ST PART:
    # Generate 1000 conformers ensembles for each ligand:
    # Get rdkit mols from the smiles under a specific uniprot id:
    mols = Chem.SmilesMolSupplier(PATH_TO_SMILES_TO_GENERATE_CONFORMER + '/' + uniprot_id + '_smiles_to_genconformer.smi', titleLine=0)
    # Open sdf writer to store the lowest energy conformer of each molecule:
    #sdf = Chem.SDWriter(PATH_TO_CONFORMERS_TO_DOCK + '/' + uniprot_id + '_minimum_energy_conformers.sdf')
    # Dictionary where all conformers are stored:
    conformers_dictionary = {}
    # List to include mols that rdkit cannot generate the conformer for
    problematic_mols = []
    t1 = time.time()
    for mol in mols:
        try:
            # Add hydrogens to molecule
            mol_Hs = AllChem.AddHs(mol)
            # Embed 1000 conformers
            list_of_ids = AllChem.EmbedMultipleConfs(mol=mol_Hs, numConfs=200, maxAttempts=200, numThreads=0,
                                                     useExpTorsionAnglePrefs=True, useBasicKnowledge=True,
                                                     useSmallRingTorsions=True, useMacrocycleTorsions=True)
            # Optimize conformers with UFF
            list_of_tuples = AllChem.UFFOptimizeMoleculeConfs(mol_Hs, maxIters=10000, numThreads=0)
            # Get the id of the lowest energy conformer
            #energy_list = []
            #for index, entry in enumerate(list_of_tuples):
            #    if entry[0] == 0:  # minimization converged
            #        energy_list.append([index, entry[1]])
            #energy_list.sort(key=lambda x: x[1])
            #id = energy_list[0][0]
            # Save the conformer
            #sdf.write(Chem.MolFromMolBlock(Chem.MolToMolBlock(mol_Hs, confId=id)))
            # Add conformers to the conformers_dictionary
            conformers_dictionary[mol.GetProp("_Name")] = mol_Hs
        except Exception as e:
            log.write('Problematic mol is ' + mol.GetProp('_Name') + '\n')
            log.write('Exception is ' + str(e) + '\n')
            #print('Problematic mol is ' + mol.GetProp('_Name'))
            #print('Exception is ' + str(e))
            problematic_mols += [mol.GetProp('_Name')]
            problem = open(PATH_TO_CONFORMERS_TO_DOCK + '/' + uniprot_id + '/' + problematic_mols[-1] + '_did_not_work.txt', 'w')
            problem.close()

    t2 = time.time()
    # How much time did the first part take?
    log.write('time_taken_to_generate_conformers: ' + str((t2 - t1) / 60) + ' minutes' + '\n' + '' + '\n')
    #print('time_taken_to_generate_conformers:', (t2-t1)/60, 'minutes')

    # 2ND PART:
    # Now from the pool of conformers choose the conformer with the least rmsd to the templates

    # Create folder where the selected conformers are saved
    os.makedirs(PATH_TO_CONFORMERS_TO_DOCK + '/' + uniprot_id, exist_ok=True)

    # Get all protein - ligand pairs corresponding to a uniprot_id
    proteins_ligands = [(key.split('_')[0], key.split('_')[1]) for key in conformers_dictionary.keys()]

    processes = []
    for pdb_id, ligand_id in proteins_ligands:
        if ligand_id in problematic_mols:
            continue
        processes += [multiprocessing.Process(target=choose_conformers, args=(uniprot_id, pdb_id, ligand_id, conformers_dictionary, problematic_mols))]

    # Run processes
    for p in processes:
        p.start()

    # Exit the completed processes
    for p in processes:
        p.join()

    tt = time.time()
    # How much time did the second part take?
    log.write('time_taken_to_run_everything: ' + str((tt - t) / 60) + ' minutes')
    log.close()
    #print('time_taken_run_everything:', (tt-t)/60, 'minutes')

def choose_conformers(uniprot_id, pdb_id, ligand_id, conformers_dictionary, problematic_mols):

    # Get model_mol
    model_mol = Chem.SDMolSupplier(PATH_TO_REFERENCE_LIGANDS_FOLDER + '/' + uniprot_id + '/' + pdb_id+'_'+ligand_id + '.sdf')[0]
    model_mol_2d = Chem.MolFromSmiles(Chem.MolToSmiles(model_mol))

    # Open sdf writer to write selected conformers
    writer = Chem.SDWriter(
        PATH_TO_CONFORMERS_TO_DOCK + '/' + uniprot_id + '/' + pdb_id + '_' + ligand_id + '_minimum_rmsd_conformers.sdf')
    # Get mols to calculate mcs
    mols = Chem.SmilesMolSupplier(PATH_TO_SMILES_TO_GENERATE_CONFORMER + '/' + uniprot_id + '_smiles_to_genconformer.smi', titleLine=0)

    for mol in mols:
        if mol.GetProp('_Name') in problematic_mols:
            continue
        try:
            # define maximum common substructure from soft constraint:
            mcs_soft = rdFMCS.FindMCS([model_mol_2d, mol], matchValences=True,
                                      ringMatchesRingOnly=True)
            mcs_strict = rdFMCS.FindMCS([model_mol_2d, mol], matchValences=True,
                                      ringMatchesRingOnly=True, completeRingsOnly=True,
                                      bondCompare=rdkit.Chem.rdFMCS.BondCompare.CompareOrderExact)
            # If FindMCS failed or if the mcs found is too small skip the molecule
            if mcs_strict.canceled == True or mcs_strict.numAtoms < 5:
                #print('No MCS shared:', mol.GetProp("_Name"))
                #print('')
                continue
            #else:
            #    print('Enough MCS shared:', mol.GetProp("_Name"))
            #    print('')
            # Get rdkit mol from the smarts strings corresponding to the strict mcs
            mcs_strict_mol = Chem.MolFromSmarts(mcs_strict.smartsString)
            # Calculate properties associated with the Smarts molecule
            mcs_strict_mol.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(mcs_strict_mol)
            # Define core
            core_strict_mol = get_core_constraint(ref_mol=model_mol, mcs_mol=mcs_strict_mol)
            # Get conformers corresponding to mol from conformers_dictionary
            mol_Hs = conformers_dictionary[mol.GetProp("_Name")]
            # Get match between conformers and the core
            match = mol_Hs.GetSubstructMatch(core_strict_mol)
            # Mapping between the atoms of the template and the conformers
            algMap = [(j, i) for i, j in enumerate(match)]
            # List with the indexes of all conformers
            lista = [conf.GetId() for conf in mol_Hs.GetConformers()]
            # List with the all conformers
            conformers = [Chem.MolFromMolBlock(Chem.MolToMolBlock(mol_Hs, confId=id)) for id in lista]
            # List with the rmsd between the mcs of template and conformers
            rmsd_s = [AllChem.AlignMol(conformers[index], core_strict_mol, atomMap=algMap) for index in range(len(conformers))]
            # Select all conformers with less than 0.5 A rmsd - this is the FlexX template based threshold
            indexes_below_threshold = [index for index in range(len(conformers)) if rmsd_s[index] < 0.5]
            if len(indexes_below_threshold) == 0:
                indexes_below_threshold = [index for index in range(len(conformers)) if np.round(rmsd_s[index],1) == np.round(min(rmsd_s),1)]
                #print('Strict mcs NOT found below 0.5 A ', min(rmsd_s), len(indexes_below_threshold))
                #conformer = conformers[rmsd_s.index(min(rmsd_s))]
                # Add name to conformer
                #conformer.SetProp("_Name", mol.GetProp("_Name"))
                # write conformer in .sdf file
                #writer.write(conformer)
            #else:
            #    print('Strict mcs found below 0.5 A ', len(indexes_below_threshold))
            # Get rdkit mol from the smarts strings corresponding to the soft mcs
            mcs_soft_mol = Chem.MolFromSmarts(mcs_soft.smartsString)
            # Calculate properties associated with the Smarts molecule
            mcs_soft_mol.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(mcs_soft_mol)
            # Define core
            core_soft_mol = get_core_constraint(ref_mol=model_mol, mcs_mol=mcs_soft_mol)
            # Get conformers corresponding to mol from conformers_dictionary
            mol_Hs = conformers_dictionary[mol.GetProp("_Name")]
            # Get match between conformers and the core
            match = mol_Hs.GetSubstructMatch(core_soft_mol)
            # Mapping between the atoms of the template and the conformers
            algMap = [(j, i) for i, j in enumerate(match)]
            # List with the rmsd between the mcs of template and conformers
            rmsd_s = [AllChem.AlignMol(conformers[index], core_soft_mol, atomMap=algMap) for index in indexes_below_threshold]
            # Selected conformer - the conformer with the minimum rmsd
            #print('rmsd to soft mcs:', min(rmsd_s))
            conformer = conformers[rmsd_s.index(min(rmsd_s))]
            # Add name to conformer
            conformer.SetProp("_Name", mol.GetProp("_Name"))
            # write conformer in .sdf file
            writer.write(conformer)
        except:
            continue

    writer.close()

def get_core_constraint(ref_mol, mcs_mol):
    ref_match = ref_mol.GetSubstructMatch(mcs_mol)

    rwmol = Chem.RWMol(mcs_mol)
    rwconf = Chem.Conformer(rwmol.GetNumAtoms())

    matches = rwmol.GetSubstructMatch(mcs_mol)

    ref_conf = ref_mol.GetConformer()
    for i, match in enumerate(matches):
        rwconf.SetAtomPosition(match, ref_conf.GetAtomPosition(ref_match[i]))

    rwmol.AddConformer(rwconf)

    return rwmol

def create_initial_reference_dictionary(path_to_reference_ligands_folder=PATH_TO_REFERENCE_LIGANDS_FOLDER):
    '''organize dictionary based on the reference ligands'''
    reference_dictionary = {}
    uniprot_ids = os.listdir(path_to_reference_ligands_folder)
    for uniprot_id in uniprot_ids:
        reference_dictionary[uniprot_id] = []
        for references in os.listdir(path_to_reference_ligands_folder + '/' + uniprot_id):
            reference_dictionary[uniprot_id] += [references]
    return reference_dictionary

def create_final_reference_dictionary(path_to_conformers_to_dock=PATH_TO_CONFORMERS_TO_DOCK):
    '''organize dictionary based on the reference ligands'''
    reference_dictionary = {}
    uniprot_ids = os.listdir(path_to_conformers_to_dock)
    uniprot_ids = [uniprot_id for uniprot_id in uniprot_ids if len(uniprot_id.split('_'))==1]
    for uniprot_id in uniprot_ids:
        reference_dictionary[uniprot_id] = []
        for references in os.listdir(path_to_conformers_to_dock + '/' + uniprot_id):
            reference_dictionary[uniprot_id] += [references]
    return reference_dictionary

if __name__ == "__main__":

    main()
