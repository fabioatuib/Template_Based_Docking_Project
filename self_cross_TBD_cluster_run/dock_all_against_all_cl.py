
# import necessary modules
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd
from modules import Mol2Writer

PATH_TO_EXTERNAL_SOFTWARE_TOOLS = '/cluster/home/fol007/external_software_tools'
PATH_TO_FLEXX = f'{PATH_TO_EXTERNAL_SOFTWARE_TOOLS}/flexx-4.5.0-Linux-x64/flexx'
PATH_TO_HYDE = f'{PATH_TO_EXTERNAL_SOFTWARE_TOOLS}/hydescorer-1.1.0-Linux-x64/hydescorer'
PATH_TO_DOCKRMSD = f'{PATH_TO_EXTERNAL_SOFTWARE_TOOLS}/DockRMSD/DockRMSD'

PATH_TO_BASE_FOLDER = '/cluster/home/fol007/template_based_docking_15_03_21'
PATH_TO_PDB_FOLDER = f'{PATH_TO_BASE_FOLDER}/pdb_files_originals'
PATH_TO_DEFINITION_FILES_FOLDER = f'{PATH_TO_BASE_FOLDER}/definition_files'
PATH_TO_REFERENCE_LIGANDS_FOLDER = f'{PATH_TO_BASE_FOLDER}/reference_ligands'
PATH_TO_CONFORMERS_FOLDER = f'{PATH_TO_BASE_FOLDER}/conformers_to_dock'
PATH_TO_DOCKED = f'{PATH_TO_BASE_FOLDER}/docked'
PATH_TO_SCORED = f'{PATH_TO_BASE_FOLDER}/scored'
PATH_TO_FILTERED = f'{PATH_TO_BASE_FOLDER}/filtered'
PATH_TO_ALIGNED = f'{PATH_TO_BASE_FOLDER}/aligned'
PATH_TO_LOG = f'{PATH_TO_BASE_FOLDER}/logs_template_based_docking'


def main():

    os.makedirs(PATH_TO_ALIGNED + '/' + 'done', exist_ok=True)
    os.makedirs(PATH_TO_ALIGNED + '/' + 'temp', exist_ok=True)

    initial_reference_dictionary = create_initial_reference_dictionary()
    final_reference_dictionary = create_final_reference_dictionary()

    print(initial_reference_dictionary)
    print(final_reference_dictionary)

    already_done_uniprot = []
    for uniprot_id in final_reference_dictionary.keys():
        if len(initial_reference_dictionary[uniprot_id]) == len(final_reference_dictionary[uniprot_id]):
            already_done_uniprot += [uniprot_id]

    uniprot_count = [(key, len(initial_reference_dictionary[key])) for key in initial_reference_dictionary.keys()]
    uniprot_count.sort(key=lambda x: x[1])

    # Do template based docking with FlexX and then scoring and then filter and then align and then rmsd
    for uniprot_id, count in uniprot_count:
        if uniprot_id in already_done_uniprot:
            continue

        os.makedirs(PATH_TO_DOCKED + '/' + uniprot_id, exist_ok=True)
        os.makedirs(PATH_TO_SCORED + '/' + uniprot_id, exist_ok=True)
        os.makedirs(PATH_TO_FILTERED + '/' + uniprot_id, exist_ok=True)
        os.makedirs(PATH_TO_ALIGNED + '/' + 'done' + '/' + uniprot_id, exist_ok=True)
        os.makedirs(PATH_TO_ALIGNED + '/' + uniprot_id, exist_ok=True)
        os.makedirs(PATH_TO_ALIGNED + '/' + uniprot_id + '/' + 'confirmed', exist_ok=True)
        os.makedirs(PATH_TO_LOG + '/' + uniprot_id, exist_ok=True)

        for ref in initial_reference_dictionary[uniprot_id]:
            if uniprot_id in final_reference_dictionary.keys():
                if ref in final_reference_dictionary[uniprot_id]:
                    continue

            if os.path.isfile(PATH_TO_CONFORMERS_FOLDER + '/' + uniprot_id + '/' + ref.split('.')[0] + '_minimum_rmsd_conformers.sdf'):

                boolean = run_sequence(uniprot_id, ref)

                if boolean == True:
                    done = open(PATH_TO_ALIGNED + '/' + 'done' + '/' + uniprot_id + '/' + ref, 'w')
                    done.close()

def run_sequence(uniprot_id, ref):

    log = open(PATH_TO_LOG + '/' + uniprot_id + '/' + ref.split('.')[0] + '.txt', 'w')
    log.write('uniprot_id: ' + uniprot_id + '\n' + '' + '\n')

    #print(uniprot_id)
    pdb_id = ref.split('_')[0]
    #print('Dock in protein:', pdb_id)
    log.write('Dock in protein: ' + pdb_id + '\n')
    path_to_protein = PATH_TO_PDB_FOLDER + '/' + uniprot_id + '/' + pdb_id + '.pdb'
    path_to_ref_lig = PATH_TO_REFERENCE_LIGANDS_FOLDER + '/' + uniprot_id + '/' + ref
    path_to_conformers = PATH_TO_CONFORMERS_FOLDER + '/' + uniprot_id + '/' + ref.split('.')[0] + '_minimum_rmsd_conformers.sdf'
    path_to_definition_file = PATH_TO_DEFINITION_FILES_FOLDER + '/' + uniprot_id + '/' + ref.split('.')[0] + '/' + ref.split('.')[0]

    try:
        # Dock
        #print('Starting Dock:')
        log.write('Starting Dock:' + '\n')
        docked = 'docked_in_' + ref.split('.')[0] + '.sdf'
        os.system(PATH_TO_FLEXX +
                  ' -i ' + path_to_conformers +
                  ' -o ' + PATH_TO_DOCKED + '/' + uniprot_id + '/' + docked +
                  ' --docking-definition ' + path_to_definition_file + '_docking_definition.ecf' +
                  ' -t ' + path_to_ref_lig +
                  ' --max-nof-conf 3')
        #print('Finished Dock:')
        log.write('Finished Dock: ' + '\n')
        path_to_docked = PATH_TO_DOCKED + '/' + uniprot_id + '/' + docked
        path_to_scored = PATH_TO_SCORED + '/' + uniprot_id + '/' + 'scored_' + docked

        # Score
        #print('Starting Score:')
        log.write('Starting Score: ' + '\n')
        os.system(PATH_TO_HYDE +
                  ' -i ' + path_to_docked +
                  ' -o ' + path_to_scored +
                  ' --binding-site-definition ' + path_to_definition_file + '_scoring_definition.ecf')
        #print('Finished Score')
        log.write('Finished Score' + '\n')

        # Filter
        print('Starting Filter')
        log.write('Starting Filter' + '\n')
        filter(path_to_scored=path_to_scored, uniprot_id=uniprot_id, docked=docked, PATH_TO_FILTERED=PATH_TO_FILTERED)
        #print('Finished Filter')
        log.write('Finished Filter' + '\n')

        # Align
        #print('Start Align:')
        log.write('Start Align:' + '\n')
        align(uniprot_id=uniprot_id, docked=docked, ref=ref, log=log, path_to_protein=path_to_protein)
        #print('Finished Align')
        log.write('Finished Align' + '\n')

        # RMSD
        #print('Start rmsd')
        log.write('Start rmsd' + '\n')
        rmsd(uniprot_id=uniprot_id, pdb_id=pdb_id, ref=ref, log=log, PATH_TO_ALIGNED=PATH_TO_ALIGNED)
        #print('Finish rmsd')
        log.write('Finish rmsd' + '\n')

        log.close()
        return True
    except Exception as e:
        #print('Exception: ' + str(e))
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

def align(uniprot_id, docked, ref, path_to_protein, log,
          PATH_TO_FILTERED=PATH_TO_FILTERED,
          PATH_TO_ALIGNED=PATH_TO_ALIGNED,
          PATH_TO_REFERENCE_LIGANDS_FOLDER=PATH_TO_REFERENCE_LIGANDS_FOLDER):

    sdf_reader = Chem.SDMolSupplier(PATH_TO_FILTERED + '/' + uniprot_id + '/' + 'filtered_' + 'scored_' + docked)

    path_to_aligned_poses = PATH_TO_ALIGNED + '/' + uniprot_id + '/' + ref.split('.')[0]
    os.makedirs(path_to_aligned_poses, exist_ok=True)
    for docked_mol in sdf_reader:
        pose_name = docked_mol.GetProp('_Name')
        path_to_pose = PATH_TO_ALIGNED + '/' + 'temp' + '/' + pose_name + '.sdf'
        writer = Chem.SDWriter(path_to_pose)
        writer.write(docked_mol)
        writer.close()

        pdb_id_of_confirmed_pose = pose_name.split('_')[0]
        ligand_name = pose_name.split('_')[1]
        path_to_confirmed = PATH_TO_PDB_FOLDER + '/' + uniprot_id + '/' + pdb_id_of_confirmed_pose + '.pdb'
        path_to_confirmed_pose = PATH_TO_REFERENCE_LIGANDS_FOLDER + '/' + uniprot_id + '/' + pdb_id_of_confirmed_pose + '_' + ligand_name + '.sdf'

        cmd.load(path_to_protein, object='docked_in')
        cmd.load(path_to_confirmed, object='confirmed')
        cmd.load(path_to_confirmed_pose, object='confirmed_pose')
        cmd.load(path_to_pose, object='pose')

        cmd.alter(selection='pose', expression='resn="pose"')

        cmd.extract('together_docked', selection='(docked_in, pose)')

        cmd.select(name='together_docked_selected', selection='br. together_docked within 10 of /together_docked///pose')
        cmd.select(name='confirmed_selected', selection='br. confirmed within 10 of confirmed_pose')

        cmd.align('together_docked_selected', 'confirmed_selected')

        path_to_aligned_complex = PATH_TO_ALIGNED + '/' + 'temp' + '/' + pose_name + '_in_' + ref.split('.')[0] + '.pdb'
        cmd.save(filename=path_to_aligned_complex, selection='together_docked')

        cmd.reinitialize()

        file = open(path_to_aligned_complex, 'r')
        pdbblock = file.read()
        file.close()

        pdbblock = pdbblock.split('\n')
        keep = []
        for i in pdbblock:
            if ('HETATM' in i) and ('pose' in i):
                keep += [i]
        keep = '\n'.join(keep)

        try:
            mol = Chem.MolFromPDBBlock(keep, sanitize=True)
            mol = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)[0]
            mol = AllChem.AssignBondOrdersFromTemplate(refmol=docked_mol, mol=mol)
            AllChem.AssignStereochemistryFrom3D(mol)
        except Exception as e:
            #print('Exception:', e)
            #print('Problem with', pose_name)
            #print('')
            log.write('Exception: ' + str(e) + '\n')
            log.write('Problem with' + pose_name + '\n')
            log.write('\n')
            continue

        Mol2Writer.MolToMol2File(mol=mol, filename=path_to_aligned_poses + '/aligned_' + pose_name + '.mol2')

        path_to_confirmed = PATH_TO_REFERENCE_LIGANDS_FOLDER + '/' + uniprot_id
        try:
            path_to_confirmed_pose = path_to_confirmed + '/' + pose_name[:-2]

            mol_conf = Chem.SDMolSupplier(path_to_confirmed_pose + '.sdf')[0]

            Mol2Writer.MolToMol2File(mol=mol_conf,
                                     filename=PATH_TO_ALIGNED + '/' + uniprot_id + '/confirmed/' + pose_name[:-2] + '.mol2')
        except:
            #print(pose_name, 'does not have corresponding confirmed pose...')
            log.write(pose_name + ' does not have corresponding confirmed pose...' + '\n')
    os.system('rm ' + PATH_TO_ALIGNED + '/' + 'temp' + '/*')

def rmsd(uniprot_id, pdb_id, ref, log, PATH_TO_ALIGNED=PATH_TO_ALIGNED):
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
            #print('did not work')
            #print(template, conf)
            #print('')
            log.write('did not work' + '\n')
            log.write(pdb_id + ' ' + conf + '\n')
            log.write('\n')
        else:
            pose = pose[0]
            path_to_docked = PATH_TO_ALIGNED + '/' + uniprot_id + '/' + template + '/' + pose
            rmsd_df['template'] += [template]
            rmsd_df['docked'] += [conf]
            rmsd = DockRMSD(path_to_docked, path_to_confirmed_pose)
            rmsd_df['rmsd'] += [rmsd]

            #print('worked')
            #print(template, conf, rmsd)
            #print('')
            log.write('worked' + '\n')
            log.write(template + ' ' + conf + ' ' + str(rmsd) + '\n')
            log.write('\n')
    rmsd_df = pd.DataFrame(rmsd_df)
    rmsd_df.to_csv(PATH_TO_ALIGNED+'/'+uniprot_id+'/'+template+'/'+'rmsd_all_in_'+template+'.csv', index=False)

def DockRMSD(query, template):
    '''returns the RMSD between the docked pose and the confirmed pose.'''
    os.system(PATH_TO_DOCKRMSD+' '+query+' '+template+' > '+PATH_TO_ALIGNED+'/temp/output.txt')
    try:
        file = open(PATH_TO_ALIGNED+'/temp/output.txt', 'r')
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
    uniprot_ids = [uniprot_id for uniprot_id in uniprot_ids if len(uniprot_id.split('_'))==1]
    for uniprot_id in uniprot_ids:
        reference_dictionary[uniprot_id] = []
        for references in os.listdir(path_to_done + '/' + uniprot_id):
            reference_dictionary[uniprot_id] += [references]
    return reference_dictionary

if __name__ == "__main__":

    main()