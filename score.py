# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 17:48:51 2023

@author: aksha
"""
import os 
import subprocess
import time
def convert_ligand_format(ligand_, new_format): 
    """Converts a ligand file to a different file format using the Open Babel tool.
        Args:
            ligand_ (str): The path to the input ligand file.
            new_format (str): The desired output format for the ligand file.
    
        Returns:
            None
    
        Raises:
            Exception: If the input file does not exist, or if the Open Babel tool is not installed.
    
        Examples:
            To convert a ligand file from mol2 format to pdbqt format:
            >>> convert_ligand_format('./ligands/ligand1.mol2', 'pdbqt')
    """
    input_format = ligand_.split('.')[-1]
    os.system('obabel {} -O {}'.format(ligand_, ligand_.replace(input_format, new_format)))

def run_vina_scoring(receptor, lig_path): 
    """
    Runs Vina scoring on the given receptor-ligand pair and returns the binding affinity score.

    Parameters:
    -----------
    receptor: str
        Path to the receptor file in pdbqt format.
    lig_path: str
        Path to the ligand file in pdbqt format.
    
    Returns:
    --------
    float
        Vina binding affinity score for the receptor-ligand pair.

    Raises:
    -------
    Exception
        If receptor file format is not pdbqt.
    Exception
        If receptor file is not found in the provided path.
    Exception
        If ligand file is not found in the provided path.
    """
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))
        
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        print('Ligand needs to be in pdbqt format. Converting ligand format using obabel.')
        convert_ligand_format(lig_path, 'pdbqt')
        lig_path = lig_path.replace(lig_format, 'pdbqt')

    cmd = ['./executables/smina', '--receptor', receptor, '-l', lig_path, '--score_only', '--scoring', 'vina']    
    command_run = subprocess.run(cmd, capture_output=True)
    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = float([x for x in command_out if 'Affinity' in x][0].split(' ')[1])

    return command_out

def check_energy(lig_): 
    """
    Check the quality of a generated structure by computing its total energy using the Open Babel obenergy tool.

    Parameters:
        lig_ (str): the name of the ligand file in PDBQT format.

    Returns:
        total_energy (float): the computed total energy of the ligand in Kcal/mol.
    """
    # Check the quality of generated structure (some post-processing quality control):
    try: 
        ob_cmd = ['obenergy', './outputs/pose_{}.pdbqt'.format(lig_.split('.')[0])]
        command_obabel_check = subprocess.run(ob_cmd, capture_output=True)
        command_obabel_check = command_obabel_check.stdout.decode("utf-8").split('\n')[-2]
        total_energy         = float(command_obabel_check.split(' ')[-2])
    except: 
        total_energy = 10000 # Calculation has failed. 
        
    return total_energy

if __name__ == '__main__':
    trial_number = 3
    # Run DiffDock 
    os.system(f'python -m inference --protein_ligand_csv data/trials/trial_{trial_number}.csv --out_dir results/trial_{trial_number} --inference_steps 20 --samples_per_complex 20 --batch_size 10 --actual_steps 18 --no_final_step_noise')
    # Change protein name as needed
    os.system('obabel ./data/1a0q/trim58.pdb -O ./data/1a0q/trim58.pdbqt')
    # Score the results of DiffDOck (Note: more than one)
    diffdock_results = os.listdir(f'./results/trial_{trial_number}')
    
    vina_scores_all = {}
    for dir_ in diffdock_results: 
        path_= f'./results/trial_{trial_number}/{dir_}/rank1.sdf' # filepath of the ligand that we want to do complex prediction on
        convert_ligand_format(path_, 'pdbqt')
        path_ = f'./results/trial_{trial_number}/{dir_}/rank1.pdbqt'
        try:
            result = check_energy(path_)
            vina_score = run_vina_scoring(receptor='./data/1a0q/trim58.pdbqt', lig_path=path_)
            vina_scores_all[dir_] = vina_score
        except Exception as e:
            vina_scores_all[dir_] = "Error in calculation"
    print(vina_scores_all)
    