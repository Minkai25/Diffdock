import os 
import subprocess
import re
import shutil

def confidence_score(filename):
    # Define a regular expression pattern to match the score in the filename.
    pattern = r"rank1_confidence([+-]?\d+\.\d{2})\.sdf"

    # Use re.search to find the match in the filename.
    match = re.search(pattern, filename)

    if match:
        # Extract the score from the matched group and convert it to a float.
        score = float(match.group(1))
        return score
    else:
        raise ValueError("Filename format is incorrect or score not found in the filename.")

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

def find_ligand_center(pdbqt_filename):
    x_sum, y_sum, z_sum = 0.0, 0.0, 0.0
    num_atoms = 0

    with open(pdbqt_filename, 'r') as pdbqt_file:
        for line in pdbqt_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Assuming the coordinates start at column 30 and are 8 characters wide for each (X, Y, Z).
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                x_sum += x
                y_sum += y
                z_sum += z
                num_atoms += 1

    if num_atoms == 0:
        raise ValueError("No ligand atoms found in the PDBQT file.")

    center_x = x_sum / num_atoms
    center_y = y_sum / num_atoms
    center_z = z_sum / num_atoms

    return round(center_x, 3), round(center_y, 3), round(center_z, 3)

if __name__ == '__main__':
    trial_number = 1
    # Run DiffDock 
    os.system(f'python -m inference --protein_ligand_csv data/trials/trial_{trial_number}.csv --out_dir results/trial_{trial_number} --inference_steps 20 --samples_per_complex 5 --batch_size 10 --actual_steps 18 --no_final_step_noise')
    # Change protein name as needed
    os.system('obabel ./data/1a0q/trim58.pdb -O ./data/1a0q/trim58.pdbqt')
    # Score the results of DiffDOck (Note: more than one)
    diffdock_results = [f for f in os.listdir(f'./results/trial_{trial_number}') if not f.startswith('.')]
    vina_scores_all = {}
    for dir_ in diffdock_results:
        path_ = None
        files = [f for f in os.listdir(f'./results/trial_{trial_number}/{dir_}') if not f.startswith('.')]
        for file in files:
            if re.match('rank1_confidence', file):
                confidence_score = file[-8:-4]
                path_ = f'./results/trial_{trial_number}/{dir_}-{confidence_score}.pdbqt'
                os.system(f'obabel ./results/trial_{trial_number}/{dir_}/{file} -O {path_}')
        shutil.rmtree(f'./results/trial_{trial_number}/{dir_}')   
        try:
            result = check_energy(path_)
            vina_score = run_vina_scoring(receptor='./data/1a0q/trim58.pdbqt', lig_path=path_)
            vina_scores_all[dir_] = vina_score
        except Exception as e:
            vina_scores_all[dir_] = "Error in calculation"
    with open(f'./results/trial_{trial_number}/dict.txt', 'w') as f:
        f.write('dict = ' + str(vina_scores_all) + '\n')  