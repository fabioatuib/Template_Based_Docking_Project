import os

PATH_TO_BASE_FOLDER = '/home/fol007/PhD_Project/Template_Based_Docking_Project/cross_free_docking_cluster_run'
PATH_TO_REFERENCE_LIGANDS_FOLDER = f'{PATH_TO_BASE_FOLDER}/reference_ligands'
PATH_TO_JOB_SCRIPTS = f'{PATH_TO_BASE_FOLDER}/job_scripts'

cpus_per_task_max = 40
correction_factor = 2
time_per_dock = 40

# Values to calculate:
cpus_per_task = None
time = None


def create_initial_reference_dictionary(path_to_reference_ligands_folder=PATH_TO_REFERENCE_LIGANDS_FOLDER):
    '''organize dictionary based on the reference ligands'''
    reference_dictionary = {}
    uniprot_ids = os.listdir(path_to_reference_ligands_folder)
    for uniprot_id in uniprot_ids:
        reference_dictionary[uniprot_id] = []
        for references in os.listdir(path_to_reference_ligands_folder + '/' + uniprot_id):
            reference_dictionary[uniprot_id] += [references.split('.')[0]]
    return reference_dictionary


def create_job_script(uniprot_id, time, cpus_per_task, size):
    job = open(PATH_TO_JOB_SCRIPTS + '/' + uniprot_id + '_' + str(size) + '.job', 'w')
    job.write('#!/bin/bash')
    job.write('\n\n')
    job.write('#######')
    job.write('\n')
    job.write('# job #')
    job.write('\n')
    job.write('#######')
    job.write('\n\n')
    job.write('## Substitute with your project name:')
    job.write('\n')
    job.write('#SBATCH --account=NN9376K')
    job.write('\n\n')
    job.write('#SBATCH --job-name=docking_' + uniprot_id)
    job.write('\n\n')
    job.write('#              d-hh:mm:ss')
    job.write('\n')
    days, hours, minutes, seconds = get_days_hours_minutes_seconds_from_seconds(time)
    job.write('#SBATCH --time=' + '%d-%02d:%02d:%02d' % (days, hours, minutes, seconds))
    job.write('\n\n')
    job.write('#SBATCH --output=slurm_outputs/docking_' + uniprot_id + '_slurm_output.out')
    job.write('\n\n')
    job.write('#SBATCH --cpus-per-task=' + str(cpus_per_task))
    job.write('\n\n')
    job.write('#SBATCH --mem-per-cpu=2G')
    job.write('\n\n')
    job.write('# you must not place bash commands before the last #SBATCH directive')
    job.write('\n\n')
    job.write('## Set safer defaults for bash')
    job.write('\n\n')
    job.write('set -o errexit')
    job.write('\n')
    job.write('set -o nounset')
    job.write('\n\n')
    job.write('module --quiet purge  # Clear any inherited modules')
    job.write('\n\n')
    job.write('python dock_all_against_all_cl.py ' + uniprot_id)
    job.close()


def get_days_hours_minutes_seconds_from_seconds(time):
    days = time // (24 * 3600)
    time = time % (24 * 3600)
    hours = time // 3600
    time %= 3600
    minutes = time // 60
    time %= 60
    seconds = time
    return days, hours, minutes, seconds


dict = create_initial_reference_dictionary()

for target in dict:
    number = len(dict[target])
    if number > cpus_per_task_max:
        cpus_per_task = cpus_per_task_max
    else:
        cpus_per_task = number
    time = (time_per_dock * number * number / cpus_per_task) * correction_factor
    create_job_script(uniprot_id=target, time=time, cpus_per_task=cpus_per_task, size=number)
