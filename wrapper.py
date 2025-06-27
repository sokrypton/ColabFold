"""
Luke Cirne
ColabFold Wrapper
Template Cycling with Experimental Distance Restraint Data
Ma Lab
"""
import os
import subprocess
import json
import getpass

def initialize_project(jobs) -> str:
    key_values = []
    # Obtain username
    try:
        username = getpass.getuser()
        print(f">>> LOGGING AS: {username}")
    except OSError:
        print(">>> COULD NOT OBTAIN USERNAME")
        print(">>> LOGGING AS: DEFAULT")
        username = "DEFUALT"
    key_values.append(("user", username))

    # Obtain job ID
    try:
        with open(jobs, "r") as file:
            print(">>> READING JSON")
            job_dict = json.load(file)
            user_JIDs = [int(job["JID"]) for job in job_dict["jobs"] if job["user"] == username]
            last_user_JID = max(user_JIDs) if user_JIDs else 0
            current_JID = last_user_JID + 1
            key_values.append(("JID", current_JID))
    except FileNotFoundError:
        current_JID = 0
        key_values.append(("JID", current_JID))

    # Obtain input file
    while True:
        input_file = input("Input file name in the format 'name.fasta': ")
        if os.path.isfile(f"./{input_file}"):
            break
        else:
            print("###### INVALID FILEPATH ######")
    key_values.append(("input_file", input_file))

    # Obtain template directory
    while True:
        temp_dir = input("Input template directory: ")
        if os.path.isdir(temp_dir):
            break
        else:
            print("###### INVALID DIRECTORY ######")
    key_values.append(("temp_dir", temp_dir))

    # Obtain num recycles
    while True:
        num_c = input("Desired number of recycles (integer) (max 5, min 0): ")
        try:
            if 0 <= int(num_c) < 6:
                break
            else:
                print("###### Invalid input ######")
        except ValueError:
                print("###### Invalid input ######")
    key_values.append(("num_recycles", num_c))    
    

    # Obtain num seeds
    while True:
        num_s = input("Desired number of seeds (integer) (min 1): ")
        try:
            if 0 < int(num_s):
                break
            else:
                print("###### Invalid input #######")
        except ValueError:
                print("###### Invalid input #######")
    key_values.append(("num_s", num_s))

    # Variable for max number of msa's
    m_e_msa = 32
    m_msa = m_e_msa // 2
    key_values.append(("m_e_msa", m_e_msa))
    key_values.append(("m_msa", m_msa))

    # Variable for full output directory by user, job id, and max msa's
    outputdir = f"{username}{current_JID}mm{m_msa}"
    key_values.append(("outputdir", outputdir))

    # Create shell script to run colabfold_batch
    current_dir = os.path.dirname(os.path.abspath(__file__))
    script_name = "wrapper.sh"
    script_path = os.path.join(current_dir, script_name)
    # Writing script
    script_content = f"""#!/bin/bash
JID={current_JID}
num_c={num_c}
seed=1
num_s={num_s}
m_e_msa={m_e_msa}
m_msa={m_msa}
inputfile=./{input_file}
outputdir={outputdir}
temp_dir={temp_dir}

colabfold_batch --pair-mode unpaired_paired --templates \\
--msa-mode mmseqs2_uniref_env \\
--custom-template-path $temp_dir \\
--max-msa $m_msa:$m_e_msa \\
--use-dropout \\
--random-seed $seed \\
--num-seeds $num_s \\
--num-recycle $num_c \\
--relax-max-iterations 100 \\
--use-gpu-relax \\
$inputfile $outputdir
    """
    # Create shell script to execute ColabFold
    with open(script_path, 'w') as file:
        print(">>> WRITING SHELL SCRIPT")
        file.write(script_content)
    append_json("jobs.json", key_values)
    return script_path


def append_json(jobs, key_values):
    # Convert list of tuples into a proper dictionary
    new_entry = {key: value for key, value in key_values}
    try:
        with open(jobs, "r") as f:
            job_dict = json.load(f)
    except FileNotFoundError:
        # Create new JSON structure if file doesn't exist
        job_dict = {"jobs": []}

    # Append the new dictionary entry
    job_dict["jobs"].append(new_entry)

    with open(jobs, "w") as f:
        print(">>> APPENDING JSON")
        json.dump(job_dict, f, indent=4)
        print("###### COMPLETE ######")


def run_colabfold(script_path):
    for run_number in range(3):
        os.chmod(script_path, 0o755)
        subprocess.run([script_path], check=True)
        filter_output(run_number)
    return 0


def filter_output(run_number):
    # Collect pdb files from output folder listed in 
    # json or directly from initialize_project()
    current_dir = os.getcwd()
    
    # loop through files and run DistanceFinder.py on each

    # sort by shortest distance and populate a temporary directory with 2 best templates
    
    return 0


def update_temp_dir(script_path, dir_name):
    with open(script_path, 'r') as file:
        lines = file.readlines()
    
    with open(script_path, 'w') as file:
        for line in lines:
            if line.startswith("temp_dir="):
                file.write(f"temp_dir={dir_name} \\")
            else:
                file.write(line)


def main():
    """ Welcome """
    print("*" * 31)
    print()
    print((" " * 7) + "ColabFold Wrapper")
    print()
    print("*" * 31)
    
    script_path = initialize_project("jobs.json")

    print(">>> ATTEMPTING TO RUN COLABFOLD")
    run_colabfold(script_path)

if __name__ == '__main__':
    main()