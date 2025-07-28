"""
Luke Cirne
ColabFold Wrapper
Template Cycling with Experimental Distance Restraint Data
Ma Lab
"""
import os
import sys
import subprocess
import shutil
import json
import getpass

def initialize_project(jobs):
    """
    Initializes a ColabFold project by gathering user input,
    creating necessary variables, generating a shell script,
    and appending metadata to a JSON log.

    Args:
        jobs (str): Path to the JSON file tracking job metadata.

    Returns:
        str: Path to the generated shell script for running ColabFold.
    """
    key_values = []
    # Obtain username
    try:
        username = getpass.getuser()
        print(f">>> LOGGING AS: {username}")
    except OSError:
        print(">>> COULD NOT OBTAIN USERNAME")
        print(">>> LOGGING AS: DEFAULT")
        username = "DEFAULT"
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
$inputfile $outputdir
    """
    # Create shell script to execute ColabFold
    with open(script_path, 'w') as file:
        print(">>> WRITING SHELL SCRIPT")
        file.write(script_content)
    append_json("jobs.json", key_values)
    return script_path


def delete_directory(dir_path):
    """
    Deletes a directory and all its contents.

    Args:
        dir_path (str): The path to the directory to delete.

    Returns:
        bool: True if the directory was successfully deleted, False otherwise.
    """
    try:
        shutil.rmtree(dir_path)
        return True
    except OSError:
        return False        


def clear_directory(dir_path):
    """
    Clears the contents of a directory without deleting the directory itself.

    Args:
        dir_path (str): Path to the directory to be cleared.
    """
    for item in os.listdir(dir_path):
        item_path = os.path.join(dir_path, item)
        if os.path.isfile(item_path):
            os.remove(item_path)  # Remove files
        elif os.path.isdir(item_path):
            shutil.rmtree(item_path)


def append_json(jobs, key_values):
    """
    Appends job metadata to a JSON file. If the file does not exist,
    it creates a new one.

    Args:
        jobs (str): Path to the JSON file.
        key_values (list): List of (key, value) tuples to append.
    """
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


def retrieve_from_current_job(jobs_file, items):
    try:
        with open(jobs_file, "r") as file:
            jobs_dict = json.load(file)
            current_job_info = jobs_dict["jobs"][-1]
            for item in items:
                match item:
                    case "user":
                        return 0
                    # TODO: finish cases and return values as an array

    except FileNotFoundError:
        print("###### JSON FILE NOT FOUND ######")
    return 0


def run_colabfold(script_path, jobs):
    """
    Runs the ColabFold shell script multiple times and filters output
    based on template distance criteria.

    Args:
        script_path (str): Path to the shell script.
        jobs (str): Path to the job metadata JSON file.
    """
    delete_directory("./recycles/")
    for run_number in range(3):
        """
        Start with three iterations for testing
        Once running, continue iterating until an ideal structure is output
        """
        os.chmod(script_path, 0o755)
        subprocess.run([script_path], check=True)
        filter_output(run_number, jobs, script_path)


def filter_output(run_number, jobs, script_path):
    """
    Filters PDB output files based on proximity to target distances.
    Updates template directory for next ColabFold iteration accordingly.

    Args:
        run_number (int): The current recycle iteration number.
        jobs (str): Path to the job metadata JSON file.
        script_path (str): Path to the ColabFold execution script.
    """
    # Load json and obtain outputdir name
    try:
        with open(jobs, "r") as f:
            jobs_dict = json.load(f)
            current_job = jobs_dict["jobs"][-1]
            current_job_values = list(current_job.values())
            outputdir = current_job_values[-1]
            temp_dir = current_job_values[3]

    except FileNotFoundError:
        print("###### JSON FILE NOT FOUND ######")
        exit()
    
    current_dir = os.getcwd()
    colabfold_output = os.listdir(f"{current_dir}/{outputdir}")

    """
    Filter to include +- 10 angstrom
    If no templates fall within range increase range to 20
    Cap at 20
    Print number of templates within range and their distances
    """
    # Loop through files and run DistanceFinder.py on each
    distances = {}
    for file in colabfold_output:
        if file.endswith(".pdb"):
            #distances.append((run_distance_finder(f"{outputdir}/{file}", "100", "473"), file))
            distances[run_distance_finder(f"{outputdir}/{file}", "100", "473")] = file

    # Filter files to only include probe distances +- 10 angstrom
    # TODO: Only include maximum of 9 files, sort by distance and take 9 best
    included_distances = {}
    for distance, filename in distances.items():
        if abs(float(distance) - 68.4) < 10:
            included_distances[distance] = filename
    # Check for any files with a distance below 10
    if not included_distances:
        for distance, filename in distances.items():
            if abs(float(distance) - 68.4) < 20:
                included_distances[distance] = filename

    # If included_distances dictionary is still empty after both checks,
    # proceed to next iteration with user provided templates 
    if not included_distances:
        print("###### NO VALID TEMPLATES PRODUCED ######")
        update_temp_dir(script_path, temp_dir)
    else:
        # Sort included_distances by shortest distance
        included_distances = dict(sorted(included_distances.items()))

        temp_dir = f"recycle{run_number + 1}"
        try:
            os.mkdir("recycles")
            print(">>> CREATING DIRECTORY FOR BEST TEMPLATES")
        except FileExistsError:
            print(">>> APPENDING TO RECYCLES DIRECTORY")

        try:
            os.mkdir(f"recycles/{temp_dir}")
        except FileExistsError:
            shutil.rmtree(f"recycles/{temp_dir}")
            os.mkdir(f"recycles/{temp_dir}")

        template_number = 0
        for distance, filename in included_distances.items():
            if template_number < 10:
                subprocess.run(["cp", f"{outputdir}/{filename}", f"recycles/{temp_dir}"])
                subprocess.run(["mv", f"recycles/{temp_dir}/{filename}", f"recycles/{temp_dir}/tmp{template_number}.pdb"])
                print(f"###### {filename} ADDED TO {temp_dir} ({distance} A) ######")
                template_number = template_number + 1

        update_temp_dir(script_path, f"recycles/{temp_dir}")
        # Clear ouput directory
        if run_number < 2:
            clear_directory(outputdir)
        #subprocess.run(["rm", "-r", outputdir])
    

def run_distance_finder(structure_file, p1, p2):
    """
    Runs an external script to calculate the distance between two residues
    in a protein structure.

    Args:
        structure_file (str): Path to the .pdb file.
        p1 (str): Residue index 1.
        p2 (str): Residue index 2.

    Returns:
        str or None: Distance in angstroms as a string, or None if failed.
    """
    distance = subprocess.run(
        ["python3", "distance_finder/DistanceFinder.py", structure_file, p1, p2],
        capture_output=True,
        text=True
    )
    if distance.returncode != 0:
        print("Error:", distance.stderr)
        return None
    return distance.stdout.strip()


def update_temp_dir(script_path, dir_name):
    """
    Updates the 'temp_dir' line in the shell script to point to a new directory.

    Args:
        script_path (str): Path to the shell script.
        dir_name (str): Name of the new template directory.
    """
    with open(script_path, 'r') as file:
        lines = file.readlines()
    
    with open(script_path, 'w') as file:
        for line in lines:
            if line.startswith("temp_dir="):
                file.write(f"temp_dir={dir_name}\n")
            else:
                file.write(line)


def main():
    """
    Entry point for the ColabFold Wrapper.
    Initializes project setup and executes the template filtering loop.
    """

    # Welcome message
    print("*" * 31)
    print()
    print((" " * 7) + "ColabFold Wrapper")
    print()
    print("*" * 31)

    jobs = "jobs.json"
    
    script_path = initialize_project(jobs)

    print(">>> ATTEMPTING TO RUN COLABFOLD\n")

    delete_directory("./recycles/")
    for run_number in range(3):
        """
        Start with three iterations for testing
        Once running, continue iterating until an ideal structure is output
        """
        os.chmod(script_path, 0o755)
        subprocess.run([script_path], check=True)
        filter_output(run_number, jobs, script_path)
    # Move recycles directory into output directory to save results
    try:
        with open(jobs, "r") as f:
            jobs_dict = json.load(f)
            current_job = jobs_dict["jobs"][-1]
            current_job_values = list(current_job.values())
            outputdir = current_job_values[-1]
    except FileNotFoundError:
        print("###### JSON FILE NOT FOUND ######")
        exit()
    subprocess.run(["mv", "./recycles/", outputdir])


if __name__ == '__main__':
    main()