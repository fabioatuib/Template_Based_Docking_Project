
import os
import time
import pandas as pd
import multiprocessing

t1 = time.time()

# Zeroth: define constants
PATH_TO_BASE_FOLDER = '/cluster/home/fol007/optimize_confirmed_with_HYDE/07_03_2021'
PATH_TO_PDB_FOLDER = f'{PATH_TO_BASE_FOLDER}/pdb_files_edited'
PATH_TO_REFERENCE_LIGANDS_FOLDER = '/cluster/home/fol007/template_based_docking_simplified/reference_ligands'
PATH_TO_OPTIMIZED_LIGANDS_FOLDER = f'{PATH_TO_BASE_FOLDER}/optimized_ligands'
PATH_TO_CSV_FILE = f'{PATH_TO_BASE_FOLDER}/expanded_Astex_with_smiles.csv'
PATH_TO_LOG = f'{PATH_TO_BASE_FOLDER}/logs'

'''
# Zeroth: define constants
PATH_TO_BASE_FOLDER = '/run/user/121976/gvfs/sftp:host=biomed3100111/Home/siv32/fol007/new_cross_docking_protocol'
PATH_TO_PDB_FOLDER = f'{PATH_TO_BASE_FOLDER}/data/pdb_files_edited'
PATH_TO_REFERENCE_LIGANDS_FOLDER = f'{PATH_TO_BASE_FOLDER}/data/reference_ligands'
PATH_TO_OPTIMIZED_LIGANDS_FOLDER = f'{PATH_TO_BASE_FOLDER}/hyde_optimize_templates/optimized_ligands'
PATH_TO_CSV_FILE = f'{PATH_TO_BASE_FOLDER}/data/expanded_Astex_with_smiles.csv'
PATH_TO_LOG = f'{PATH_TO_BASE_FOLDER}/hyde_optimize_templates/logs'
'''

# First: start by opening the .csv file that serves as a guide
df = pd.read_csv(PATH_TO_CSV_FILE)
df = df[['Uniprot_ID', 'Protein_ID', 'Ligand_Name']]
df = df.iloc[0:10]

os.makedirs(PATH_TO_LOG, exist_ok=True)
os.makedirs(PATH_TO_OPTIMIZED_LIGANDS_FOLDER, exist_ok=True)

def HYDE(path_to_reference_ligand, path_to_optimized_ligand, path_to_pdb_file, path_to_log_file):
    #os.system('/cluster/home/fol007/BioSolveIT_latest_releases/hydescorer-1.0.1-Linux-x64/hydescorer' +
    os.system('/cluster/home/fol007/hydescorer-1.1.0-Linux-x64/hydescorer' +
                  ' -i ' + path_to_reference_ligand +
                  ' -o ' + path_to_optimized_ligand +
                  ' -p ' + path_to_pdb_file +
                  ' -r ' + path_to_reference_ligand +
                  ' --thread-count' + ' 1 ' + ' > ' + path_to_log_file)

# Second: iterate over the dataframe and make every HYDE optimization run in parallel with os.system("nohup ... &")
list_of_lists = df.values.tolist()
list_of_lists = [list_of_lists[30*i:30+30*i] for i in range((len(list_of_lists)//30)-1)] + [list_of_lists[30*((len(list_of_lists)//30)-1):]]

for lista in list_of_lists:
    processes = []
    for unp, protein, ligand in lista:

        if os.path.isdir(PATH_TO_OPTIMIZED_LIGANDS_FOLDER+'/'+unp) == False:
            os.makedirs(PATH_TO_OPTIMIZED_LIGANDS_FOLDER+'/'+unp, exist_ok=True)

        path_to_pdb_file = PATH_TO_PDB_FOLDER + '/' + unp + '/' + protein + '_' + ligand + '.pdb'
        path_to_reference_ligand = PATH_TO_REFERENCE_LIGANDS_FOLDER + '/' + unp + '/' + protein + '_' + ligand + '.sdf'
        path_to_optimized_ligand = PATH_TO_OPTIMIZED_LIGANDS_FOLDER + '/' + unp + '/' + protein + '_' + ligand + '.sdf'
        path_to_log_file = PATH_TO_LOG + '/' + protein + '_' + ligand + '.log'
    
        processes += [multiprocessing.Process(target=HYDE, args=(path_to_reference_ligand, path_to_optimized_ligand, path_to_pdb_file, path_to_log_file))]

    # Run processes
    for p in processes:
        p.start()

    # Exit the completed processes
    for p in processes:
        p.join()

t2 = time.time()

time_log = open(f'{PATH_TO_BASE_FOLDER}/time_log.txt', 'w')
time_log.write(str(t2-t1) + ' seconds')
time_log.close()
