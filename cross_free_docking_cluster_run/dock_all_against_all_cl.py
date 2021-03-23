# import necessary modules
import os
import sys
import time
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd
import Mol2Writer

PATH_TO_EXTERNAL_SOFTWARE_TOOLS = '/home/fol007/software_tools'
PATH_TO_FLEXX = f'{PATH_TO_EXTERNAL_SOFTWARE_TOOLS}/flexx-4.5.0-Linux-x64/flexx'
PATH_TO_HYDE = f'{PATH_TO_EXTERNAL_SOFTWARE_TOOLS}/hydescorer-1.1.0-Linux-x64/hydescorer'
PATH_TO_DOCKRMSD = f'{PATH_TO_EXTERNAL_SOFTWARE_TOOLS}/DockRMSD/DockRMSD'

PATH_TO_BASE_FOLDER = '/home/fol007/PhD_Project/Template_Based_Docking_Project/cross_free_docking_cluster_run'
PATH_TO_PDB_FOLDER = f'{PATH_TO_BASE_FOLDER}/pdb_files_edited'
PATH_TO_REFERENCE_LIGANDS_FOLDER = f'{PATH_TO_BASE_FOLDER}/reference_ligands'
PATH_TO_CONFORMERS_FOLDER = f'{PATH_TO_BASE_FOLDER}/conformers_to_dock'
PATH_TO_DOCKED = f'{PATH_TO_BASE_FOLDER}/docked'
PATH_TO_SCORED = f'{PATH_TO_BASE_FOLDER}/scored'
PATH_TO_FILTERED = f'{PATH_TO_BASE_FOLDER}/filtered'
PATH_TO_ALIGNED = f'{PATH_TO_BASE_FOLDER}/aligned'
PATH_TO_LOG = f'{PATH_TO_BASE_FOLDER}/logs_template_based_docking'


def main():
    # create folders
    os.makedirs(PATH_TO_ALIGNED + '/' + 'done', exist_ok=True)
    os.makedirs(PATH_TO_ALIGNED + '/' + 'temp', exist_ok=True)

    initial_reference_dictionary = create_initial_reference_dictionary()
    final_reference_dictionary = create_final_reference_dictionary()

    if len(sys.argv) > 1:
        uniprot_id = sys.argv[1]
        initial_reference_dictionary = {uniprot_id:initial_reference_dictionary[uniprot_id]}
        if uniprot_id in final_reference_dictionary:
            final_reference_dictionary = {uniprot_id: final_reference_dictionary[uniprot_id]}
        else:
            final_reference_dictionary = {}

    print(initial_reference_dictionary)
    print(final_reference_dictionary)

    for uniprot_id in final_reference_dictionary.keys():
        # if they are the same length they everything is done regarding that target
        if len(initial_reference_dictionary[uniprot_id]) == len(final_reference_dictionary[uniprot_id]):
            initial_reference_dictionary.pop(uniprot_id)
            continue
        # if some are missing, remove the ones that were already done
        for already_done in final_reference_dictionary[uniprot_id]:
            initial_reference_dictionary[uniprot_id].remove(already_done)

    # order from the protein target that has more entries to dock to the one with least entries
    uniprot_count = [(key, len(initial_reference_dictionary[key])) for key in initial_reference_dictionary.keys()]
    uniprot_count.sort(key=lambda x: x[1])

    # Do free docking with FlexX and then scoring with HyDE and then filter and then align and then rmsd
    for uniprot_id, count in uniprot_count:

        os.makedirs(PATH_TO_DOCKED + '/' + uniprot_id, exist_ok=True)
        os.makedirs(PATH_TO_SCORED + '/' + uniprot_id, exist_ok=True)
        os.makedirs(PATH_TO_FILTERED + '/' + uniprot_id, exist_ok=True)
        os.makedirs(PATH_TO_ALIGNED + '/' + 'done' + '/' + uniprot_id, exist_ok=True)
        os.makedirs(PATH_TO_ALIGNED + '/' + uniprot_id, exist_ok=True)
        os.makedirs(PATH_TO_ALIGNED + '/' + uniprot_id + '/' + 'confirmed', exist_ok=True)
        os.makedirs(PATH_TO_LOG + '/' + uniprot_id, exist_ok=True)

        for ref in initial_reference_dictionary[uniprot_id]:

            if os.path.isfile(PATH_TO_CONFORMERS_FOLDER + '/' + uniprot_id + '/' + ref.split('.')[
                0] + '_minimum_rmsd_conformers.sdf') and \
                    os.path.isfile(PATH_TO_PDB_FOLDER + '/' + uniprot_id + '/' + ref.split('.')[0] + '.pdb'):

                conformers = [m.GetProp('_Name') for m in
                              Chem.SDMolSupplier(PATH_TO_CONFORMERS_FOLDER + '/' + uniprot_id + '/' + ref.split('.')[
                                  0] + '_minimum_rmsd_conformers.sdf')]
                log_df = pd.DataFrame(data=[[None] * 7] * len(conformers),
                                      columns=['template', 'to_dock', 'docked', 'scored', 'align_error', 'rmsd',
                                               'general_error'])
                log_df['template'] = ref.split('.')[0]
                log_df['to_dock'] = conformers

                boolean = run_sequence(uniprot_id, ref, log_df)

                log_df.to_csv(PATH_TO_ALIGNED + '/' + 'done' + '/' + uniprot_id + '/' + ref.split('.')[0] + '.csv',
                              index=False)

            else:
                done = open(PATH_TO_ALIGNED + '/' + 'done' + '/' + uniprot_id + '/' + ref.split('.')[0] + '.csv', 'w')
                done.close()


def run_sequence(uniprot_id, ref, log_df):
    log = open(PATH_TO_LOG + '/' + uniprot_id + '/' + ref.split('.')[0] + '.txt', 'w')
    log.write('uniprot_id: ' + uniprot_id + '\n' + '' + '\n')

    print(uniprot_id)
    pdb_id = ref.split('_')[0]
    print('Dock in protein:', pdb_id)
    log.write('Dock in protein: ' + ref.split('.')[0] + '\n')
    path_to_protein = PATH_TO_PDB_FOLDER + '/' + uniprot_id + '/' + ref.split('.')[0] + '.pdb'
    path_to_ref_lig = PATH_TO_REFERENCE_LIGANDS_FOLDER + '/' + uniprot_id + '/' + ref
    path_to_conformers = PATH_TO_CONFORMERS_FOLDER + '/' + uniprot_id + '/' + ref.split('.')[
        0] + '_minimum_rmsd_conformers.sdf'

    try:
        # Dock
        print('Starting Dock:')
        log.write('Starting Dock:' + '\n')
        docked = 'docked_in_' + ref.split('.')[0] + '.sdf'
        path_to_docked = PATH_TO_DOCKED + '/' + uniprot_id + '/' + docked
        os.system(PATH_TO_FLEXX +
                  ' -i ' + path_to_conformers +
                  ' -o ' + path_to_docked +
                  ' -r ' + path_to_ref_lig +
                  ' -p ' + path_to_protein +
                  ' --max-nof-conf 3')
        print('Finished Dock:')
        log.write('Finished Dock: ' + '\n')

        # write down result on log_df
        results = {m.GetProp('_Name')[:-2] for m in
                   Chem.SDMolSupplier(path_to_docked)}
        log_df['docked'] = 0
        for result in results:
            log_df.loc[log_df['to_dock'] == result, 'docked'] = 1

        # Score
        print('Starting Score:')
        log.write('Starting Score: ' + '\n')
        path_to_scored = PATH_TO_SCORED + '/' + uniprot_id + '/' + 'scored_' + docked
        os.system(PATH_TO_HYDE +
                  ' -i ' + path_to_docked +
                  ' -o ' + path_to_scored +
                  ' -r ' + path_to_ref_lig +
                  ' -p ' + path_to_protein)
        print('Finished Score')
        log.write('Finished Score' + '\n')

        # write down result on log_df
        results = {m.GetProp('_Name')[:-2] for m in
                   Chem.SDMolSupplier(path_to_scored)}
        log_df['scored'] = 0
        for result in results:
            log_df.loc[log_df['to_dock'] == result, 'scored'] = 1

        # Filter
        print('Starting Filter')
        log.write('Starting Filter' + '\n')
        filter(path_to_scored=path_to_scored, uniprot_id=uniprot_id, docked=docked, PATH_TO_FILTERED=PATH_TO_FILTERED)
        print('Finished Filter')
        log.write('Finished Filter' + '\n')

        # Align
        print('Start Align:')
        log.write('Start Align:' + '\n')
        align(uniprot_id=uniprot_id, docked=docked, ref=ref, log=log, log_df=log_df, path_to_protein=path_to_protein)
        print('Finished Align')
        log.write('Finished Align' + '\n')

        # RMSD
        print('Start rmsd')
        log.write('Start rmsd' + '\n')
        rmsd(uniprot_id=uniprot_id, pdb_id=pdb_id, ref=ref, log=log, log_df=log_df, PATH_TO_ALIGNED=PATH_TO_ALIGNED)
        print('Finish rmsd')
        log.write('Finish rmsd' + '\n')

        log.close()
        return True
    except Exception as e:
        log_df['general_error'] = str(e)
        print('Exception: ' + str(e))
        log.write('Exception: ' + str(e) + '\n')
        log.close()
        return False


def filter(path_to_scored, uniprot_id, docked, PATH_TO_FILTERED=PATH_TO_FILTERED):
    '''filter the scored compounds: only keep the one with the best Hyde score'''
    mols = Chem.SDMolSupplier(path_to_scored)
    mols_iter = iter([mol for mol in mols])

    mols_list = []

    boolean = True
    mol = next(mols_iter, 'STOP')
    while boolean:
        _name = mol.GetProp("_Name")
        mols_list += [[[mol.GetProp("_Name"), mol]]]

        mol = next(mols_iter, 'STOP')
        if mol == 'STOP':
            break
        while mol.GetProp("_Name")[:8] == _name[:8]:
            mols_list[-1] += [[mol.GetProp("_Name"), mol]]
            mol = next(mols_iter, 'STOP')
            if mol == 'STOP':
                boolean = False
                break

    for poses in mols_list:
        for pose in poses:
            if 'BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]' in pose[1].GetPropsAsDict():
                pose += [pose[1].GetPropsAsDict()['BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]']]
            else:
                pose += [None]
            if 'docking-score' in pose[1].GetPropsAsDict():
                pose += [pose[1].GetPropsAsDict()['docking-score']]
            else:
                pose += [None]

    new_mols_list = []

    for poses in mols_list:
        booleans = [pose[2] == None for pose in poses]
        if False in booleans:
            poses = [pose for pose in poses if pose[2] != None]
            poses.sort(key=lambda x: x[2])
            poses = [pose for pose in poses if pose[2] == poses[0][2]]
            if len(poses) == 1:
                new_mols_list += [poses[0]]
            else:
                poses.sort(key=lambda x: x[3])
                new_mols_list += [poses[0]]
        else:
            poses.sort(key=lambda x: x[3])
            new_mols_list += [poses[0]]

    file = Chem.SDWriter(PATH_TO_FILTERED + '/' + uniprot_id + '/' + 'filtered_' + 'scored_' + docked)
    for mol in new_mols_list:
        file.write(mol[1])
    file.close()


def align(uniprot_id, docked, ref, path_to_protein, log, log_df,
          PATH_TO_FILTERED=PATH_TO_FILTERED,
          PATH_TO_ALIGNED=PATH_TO_ALIGNED,
          PATH_TO_REFERENCE_LIGANDS_FOLDER=PATH_TO_REFERENCE_LIGANDS_FOLDER):
    # TODO refactor align function!
    sdf_reader = Chem.SDMolSupplier(PATH_TO_FILTERED + '/' + uniprot_id + '/' + 'filtered_' + 'scored_' + docked)

    path_to_aligned_poses = PATH_TO_ALIGNED + '/' + uniprot_id + '/' + ref.split('.')[0]
    os.makedirs(path_to_aligned_poses, exist_ok=True)
    for docked_mol in sdf_reader:
        try:
            pose_name = docked_mol.GetProp('_Name')
            path_to_pose = PATH_TO_ALIGNED + '/' + 'temp' + '/' + pose_name + '.sdf'
            writer = Chem.SDWriter(path_to_pose)
            writer.write(docked_mol)
            writer.close()

            pdb_id_of_confirmed_pose = '_'.join(pose_name.split('_')[:2])
            path_to_confirmed = PATH_TO_PDB_FOLDER + '/' + uniprot_id + '/' + pdb_id_of_confirmed_pose + '.pdb'
            path_to_confirmed_pose = PATH_TO_REFERENCE_LIGANDS_FOLDER + '/' + uniprot_id + '/' + pdb_id_of_confirmed_pose + '.sdf'

            cmd.reinitialize()

            cmd.load(path_to_protein, object='docked_in')
            cmd.load(path_to_confirmed, object='confirmed')
            cmd.load(path_to_confirmed_pose, object='confirmed_pose')
            cmd.load(path_to_pose, object='pose')

            cmd.alter(selection='pose', expression='resn="pose"')

            cmd.extract('together_docked', selection='(docked_in, pose)')

            cmd.select(name='together_docked_selected',
                       selection='br. together_docked within 10 of /together_docked///pose')
            cmd.select(name='confirmed_selected', selection='br. confirmed within 10 of confirmed_pose')

            cmd.align('together_docked_selected', 'confirmed_selected')

            path_to_aligned_complex = PATH_TO_ALIGNED + '/' + 'temp' + '/' + pose_name + '_in_' + ref.split('.')[
                0] + '.pdb'
            cmd.save(filename=path_to_aligned_complex, selection='together_docked')

            file = open(path_to_aligned_complex, 'r')
            pdbblock = file.read()
            file.close()

            pdbblock = pdbblock.split('\n')
            keep = []
            for i in pdbblock:
                if ('HETATM' in i) and ('pose' in i):
                    keep += [i]
            keep = '\n'.join(keep)

            mol = Chem.MolFromPDBBlock(keep, sanitize=True)
            mol = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)[0]
            mol = AllChem.AssignBondOrdersFromTemplate(refmol=docked_mol, mol=mol)
            AllChem.AssignStereochemistryFrom3D(mol)

            Mol2Writer.MolToMol2File(mol=mol, filename=path_to_aligned_poses + '/aligned_' + pose_name + '.mol2')

            if not os.path.isfile(PATH_TO_ALIGNED + '/' + uniprot_id + '/confirmed/' + pose_name[:-2] + '.mol2'):
                path_to_confirmed_pose = PATH_TO_REFERENCE_LIGANDS_FOLDER + '/' + uniprot_id + '/' + pose_name[
                                                                                                     :-2] + '.sdf'
                mol_conf = Chem.SDMolSupplier(path_to_confirmed_pose)[0]
                Mol2Writer.MolToMol2File(mol=mol_conf,
                                         filename=PATH_TO_ALIGNED + '/' + uniprot_id + '/confirmed/' + pose_name[
                                                                                                       :-2] + '.mol2')
        except Exception as e:
            print('Exception:', e)
            print('Problem with', pose_name)
            print('')
            log.write('Exception: ' + str(e) + '\n')
            log.write('Problem with' + pose_name + '\n')
            log.write('\n')
            log_df.loc[log_df['to_dock'] == pose_name[:-2], 'align_error'] = 'Exception: ' + str(e)
            continue

    os.system('rm ' + PATH_TO_ALIGNED + '/' + 'temp' + '/*')


def rmsd(uniprot_id, pdb_id, ref, log, log_df, PATH_TO_ALIGNED=PATH_TO_ALIGNED):
    '''Calculate rmsd'''

    template = ref.split('.')[0]
    rmsd_df = {'template': [], 'docked': [], 'rmsd': []}

    aligned_dictionary = {}
    aligned_dictionary[uniprot_id] = {}
    for folder in os.listdir(PATH_TO_ALIGNED + '/' + uniprot_id):
        if folder in ['confirmed', template]:
            aligned_dictionary[uniprot_id][folder] = []
            for mol in os.listdir(PATH_TO_ALIGNED + '/' + uniprot_id + '/' + folder):
                aligned_dictionary[uniprot_id][folder] += [mol]

    confirmed = aligned_dictionary[uniprot_id]['confirmed']

    for conf in confirmed:
        path_to_confirmed_pose = PATH_TO_ALIGNED + '/' + uniprot_id + '/confirmed/' + conf
        conf = conf.split('.')[0]

        pose = [pose for pose in aligned_dictionary[uniprot_id][template] if 'aligned_' + conf in pose]
        if pose == []:
            print('did not work')
            print(template, conf)
            print('')
            log.write('did not work' + '\n')
            log.write(pdb_id + ' ' + conf + '\n')
            log.write('\n')
            log_df.loc[log_df['to_dock'] == conf, 'rmsd'] = None
        else:
            pose = pose[0]
            path_to_docked = PATH_TO_ALIGNED + '/' + uniprot_id + '/' + template + '/' + pose
            rmsd_df['template'] += [template]
            rmsd_df['docked'] += [conf]
            rmsd = DockRMSD(path_to_docked, path_to_confirmed_pose)
            rmsd_df['rmsd'] += [rmsd]
            log_df.loc[log_df['to_dock'] == conf, 'rmsd'] = rmsd
            print('worked')
            print(template, conf, rmsd)
            print('')
            log.write('worked' + '\n')
            log.write(template + ' ' + conf + ' ' + str(rmsd) + '\n')
            log.write('\n')
    rmsd_df = pd.DataFrame(rmsd_df)
    rmsd_df.to_csv(PATH_TO_ALIGNED + '/' + uniprot_id + '/' + template + '/' + 'rmsd_all_in_' + template + '.csv',
                   index=False)


def DockRMSD(query, template):
    '''returns the RMSD between the docked pose and the confirmed pose.'''
    os.system(PATH_TO_DOCKRMSD + ' ' + query + ' ' + template + ' > ' + PATH_TO_ALIGNED + '/temp/output.txt')
    try:
        file = open(PATH_TO_ALIGNED + '/temp/output.txt', 'r')
        for line in file:
            if 'Calculated Docking RMSD:' in line:
                line = line.split('\n')[0]
                line = line.split(': ')[1]
                return float(line)
    except:
        return None


def create_initial_reference_dictionary(path_to_reference_ligands_folder=PATH_TO_REFERENCE_LIGANDS_FOLDER):
    '''organize dictionary based on the reference ligands'''
    reference_dictionary = {}
    uniprot_ids = os.listdir(path_to_reference_ligands_folder)
    for uniprot_id in uniprot_ids:
        reference_dictionary[uniprot_id] = []
        for references in os.listdir(path_to_reference_ligands_folder + '/' + uniprot_id):
            reference_dictionary[uniprot_id] += [references]
    return reference_dictionary


def create_final_reference_dictionary(path_to_done=PATH_TO_ALIGNED + '/' + 'done'):
    '''organize dictionary based on the reference ligands'''
    reference_dictionary = {}
    uniprot_ids = os.listdir(path_to_done)
    uniprot_ids = [uniprot_id for uniprot_id in uniprot_ids if len(uniprot_id.split('_')) == 1]
    for uniprot_id in uniprot_ids:
        reference_dictionary[uniprot_id] = []
        for references in os.listdir(path_to_done + '/' + uniprot_id):
            reference_dictionary[uniprot_id] += [references.split('.')[0]+'.sdf']
    return reference_dictionary


if __name__ == "__main__":
    time1=time.time()
    main()
    time2=time.time()
    print(time2-time1)
