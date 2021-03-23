
# needed modules
import os
import json
from zipfile import ZipFile
from biopandas.pdb import PandasPdb
from rdkit import Chem
import pandas as pd

# needed input
'''
path_to_pdb_file = '../../data/pdb_files_edited/P0A6I3/1SQ5_ADP.pdb'
complex_id = '1SQ5_ADP'
path_to_reference = '../../data/optimized_ligands/P0A6I3/1SQ5_ADP.sdf'
path_to_store_definition_file = '../dock_and_score_definitions/'+complex_id
'''

def generate_definitions(complex_id, path_to_pdb_file, path_to_reference, path_to_store_definition_file):
    # get structure
    ppdb = PandasPdb().read_pdb(path_to_pdb_file)

    atom_df = ppdb.df['ATOM']
    hetatm_df = ppdb.df['HETATM']

    reference = Chem.SDMolSupplier(path_to_reference)[0]
    reference = reference.GetConformer()

    # calculate distances between reference ligand atoms to all other atoms in the structure
    i = 0
    series_atoms = {}
    series_hetatoms = {}
    for position in reference.GetPositions():
        series_atoms[i] = ppdb.distance(xyz=position, records=('ATOM'))
        series_hetatoms[i] = ppdb.distance(xyz=position, records=('HETATM'))
        i += 1
    # concatenate into a dataframe
    series_atoms = pd.concat(series_atoms, axis=1)
    series_hetatoms = pd.concat(series_hetatoms, axis=1)
    # assign to each atom the minimum distance to the reference ligand
    series_atoms = series_atoms.min(axis=1)
    series_hetatoms = series_hetatoms.min(axis=1)
    # add column with distances to the 'ATOM' biopandas dataframe
    atom_df = pd.concat([series_atoms, atom_df], axis=1)
    # add column with distances to the 'HETATM' biopandas dataframe
    hetatm_df = pd.concat([series_hetatoms, hetatm_df], axis=1)

    # get residue number and chain_id of aa residues that are within 15 A
    residue_numbers_atom = atom_df.loc[atom_df[0]<15][[0, 'residue_number', 'chain_id', 'residue_name']].groupby(['residue_number', 'chain_id', 'residue_name']).min()
    # coordinates of C_alpha of aa that are within 10.5 A
    coordinates_C_alphas = []
    for ids, dist in residue_numbers_atom.itertuples():
        residue_number, chain_id, residue_name = ids[0], ids[1], ids[2]
        coordinates_C_alpha = atom_df.loc[(atom_df['atom_name']=='C') &
                                            (atom_df['residue_number']==residue_number) &
                                            (atom_df['chain_id']==chain_id) &
                                            (atom_df['residue_name']==residue_name)][
                                                                            ['x_coord', 'y_coord', 'z_coord']].iloc[:1].values.tolist()
        if len(coordinates_C_alpha) == 1:
            coordinates_C_alphas += coordinates_C_alpha
            coordinates_C_alphas[-1] = [dist] + coordinates_C_alphas[-1]

    coordinates_C_alphas = pd.DataFrame(coordinates_C_alphas, columns=[0, 'x_coord', 'y_coord', 'z_coord'])

    # coordinates of the first atom entries of other heterogens included in the binding site definition
    coordinates_first_hetatm_entry = hetatm_df.loc[(~hetatm_df['residue_name'].isin(['HOH'])) & (hetatm_df[0]<15)][
                                                [0, 'residue_number', 'x_coord', 'y_coord', 'z_coord']]
    # put indexes in a column
    coordinates_first_hetatm_entry = coordinates_first_hetatm_entry.reset_index()
    # get indexes of first atom entry
    indexes = coordinates_first_hetatm_entry[['residue_number', 'index']].groupby(['residue_number']).min()['index'].to_list()
    # coordinates of first atom entry
    coordinates_first_hetatm_entry = coordinates_first_hetatm_entry.loc[coordinates_first_hetatm_entry['index'].isin(indexes)][
                                                [0,'x_coord', 'y_coord', 'z_coord']]

    # create the docking definition .json file
    '''
    {'MaximumOverlapVolumeForDocking': '2.9',
     'bindingSite': {'definitionSpheres': [
         {'r': 0.1, 'x': 60.075, 'y': 44.607, 'z': 20.399}]},
     'dataType': 1,
     'extendedBindingSite': {'definitionSpheres': [
       {'r': 0.1, 'x': 64.574, 'y': 38.854, 'z': 21.33}]},
     'version': 3}
    '''

    docking_definition = {"MaximumOverlapVolumeForDocking": "2.9",
                          "bindingSite": {"definitionSpheres": []},
                          "dataType": 1,
                          "extendedBindingSite": {"definitionSpheres": []},
                          "version": 3}

    # add bindingSite definitionSpheres:
    for x, y, z in coordinates_C_alphas.loc[coordinates_C_alphas[0]<=10][['x_coord', 'y_coord', 'z_coord']].values:
        docking_definition['bindingSite']['definitionSpheres'] += [{'r': 0.1, 'x': x, 'y': y, 'z': z}]

    for x, y, z in coordinates_first_hetatm_entry.loc[coordinates_first_hetatm_entry[0]<=10][['x_coord', 'y_coord', 'z_coord']].values:
        docking_definition['bindingSite']['definitionSpheres'] += [{'r': 0.1, 'x': x, 'y': y, 'z': z}]

    # add extendedBindingSite definitionSpheres:
    for x, y, z in coordinates_C_alphas.loc[coordinates_C_alphas[0]<=15][['x_coord', 'y_coord', 'z_coord']].values:
        docking_definition['extendedBindingSite']['definitionSpheres'] += [{'r': 0.1, 'x': x, 'y': y, 'z': z}]

    for x, y, z in coordinates_first_hetatm_entry.loc[coordinates_first_hetatm_entry[0]<=15][['x_coord', 'y_coord', 'z_coord']].values:
        docking_definition['extendedBindingSite']['definitionSpheres'] += [{'r': 0.1, 'x': x, 'y': y, 'z': z}]


    # create the scoring definition .json file
    '''
    {'bindingSite': {'definitionSpheres': [
       {'r': 0.1, 'x': 53.102, 'y': 55.541, 'z': 14.728},
       ]},
     'dataType': 3,
     'extendedBindingSite': {'definitionSpheres': [
       {'r': 0.1, 'x': 47.396, 'y': 44.034, 'z': 24.396}]},
     'version': 3}
    '''

    scoring_definition = {"bindingSite": {"definitionSpheres": []},
                          "dataType": 3,
                          "extendedBindingSite": {"definitionSpheres": []},
                          "version": 3}

    # add bindingSite definitionSpheres:
    for x, y, z in coordinates_C_alphas.loc[coordinates_C_alphas[0]<=10][['x_coord', 'y_coord', 'z_coord']].values:
        scoring_definition['bindingSite']['definitionSpheres'] += [{'r': 0.1, 'x': x, 'y': y, 'z': z}]

    for x, y, z in coordinates_first_hetatm_entry.loc[coordinates_first_hetatm_entry[0]<=10][['x_coord', 'y_coord', 'z_coord']].values:
        scoring_definition['bindingSite']['definitionSpheres'] += [{'r': 0.1, 'x': x, 'y': y, 'z': z}]

    # add extendedBindingSite definitionSpheres:
    for x, y, z in coordinates_C_alphas.loc[coordinates_C_alphas[0]<=15][['x_coord', 'y_coord', 'z_coord']].values:
        scoring_definition['extendedBindingSite']['definitionSpheres'] += [{'r': 0.1, 'x': x, 'y': y, 'z': z}]

    for x, y, z in coordinates_first_hetatm_entry.loc[coordinates_first_hetatm_entry[0]<=15][['x_coord', 'y_coord', 'z_coord']].values:
        scoring_definition['extendedBindingSite']['definitionSpheres'] += [{'r': 0.1, 'x': x, 'y': y, 'z': z}]

    # save the definition files
    with open(path_to_store_definition_file+'/docking_definition.json', 'w') as outfile:
        json.dump(docking_definition, outfile, indent=4)

    with open(path_to_store_definition_file+'/scoring_definition.json', 'w') as outfile:
        json.dump(scoring_definition, outfile, indent=4)

    # save back to .pdb file
    ppdb.to_pdb(path_to_store_definition_file+'/protein.pdb', records={'ATOM', 'HETATM'})

    with ZipFile(path_to_store_definition_file+'/'+complex_id+'_docking_definition.ecf', 'w') as zip:
        zip.write(path_to_store_definition_file+'/protein.pdb', 'protein.pdb')
        zip.write(path_to_store_definition_file+'/docking_definition.json', 'definition.json')

    with ZipFile(path_to_store_definition_file+'/'+complex_id+'_scoring_definition.ecf', 'w') as zip:
        zip.write(path_to_store_definition_file+'/protein.pdb', 'protein.pdb')
        zip.write(path_to_store_definition_file+'/scoring_definition.json', 'definition.json')

    os.remove(path_to_store_definition_file+'/protein.pdb')
    os.remove(path_to_store_definition_file+'/docking_definition.json')
    os.remove(path_to_store_definition_file+'/scoring_definition.json')