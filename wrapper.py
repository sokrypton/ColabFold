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

def initialize_project(jobs) -> str:
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


def clear_directory(dir_path):
    try:
        shutil.rmtree(dir_path)
        return True
    except OSError:
        return False        


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


def run_colabfold(script_path, jobs):
    clear_directory("./recycles/")
    for run_number in range(3):
        """
        Start with three iterations for testing
        Once running, continue iterating until an ideal structure is output
        """
        os.chmod(script_path, 0o755)
        subprocess.run([script_path], check=True)
        filter_output(run_number, jobs, script_path)
    return 0


def filter_output(run_number, jobs, script_path):
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
        return -1
    
    current_dir = os.getcwd()
    colabfold_output = os.listdir(f"{current_dir}/{outputdir}")

    """
    TODO: Modify filter to include +- 10 angstrom
    If no templates fall within range increase range to 20
    Cap at 20
    Print number of templates within range and their distances
    """
    # Loop through files and run DistanceFinder.py on each
    distances = {}
    for file in colabfold_output:
        if file.endswith(".pdb"):
            #distances.append((run_distance_finder(f"{outputdir}/{file}", "100", "473"), file))
            distances[file] = run_distance_finder(f"{outputdir}/{file}", "100", "473")

    # Filter files to only include probe distances +- 10 angstrom
    included_distances = {}
    for filename, distance in distances.items():
        if abs(float(distance) - 68.4) < 10:
            included_distances[filename] = distance
    # Check for any files with a distance below 10
    if not included_distances:
        for filename, distance in distances.items():
            if abs(float(distance) - 68.4) < 20:
                included_distances[filename] = distance

    # If included_distances dictionary is still empty after both checks,
    # proceed to next iteration with user provided templates 
    if not included_distances:
        print("###### NO VALID TEMPLATES PRODUCED ######")
        update_temp_dir(script_path, temp_dir)
    else:
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

        template_number = 1
        for filename, distance in included_distances.items():
            subprocess.run(["cp", f"{outputdir}/{filename}", f"recycles/{temp_dir}"])
            subprocess.run(["mv", f"recycles/{temp_dir}/{filename}", f"recycles/{temp_dir}/tmp{template_number}.pdb"])
            print(f"###### {filename} ADDED TO {temp_dir} ({distance} A) ######")
            template_number = template_number + 1

            
            
        update_temp_dir(script_path, f"recycles/{temp_dir}")
        # Clear ouput directory
        #subprocess.run(["rm", "-r", outputdir])
    

def run_distance_finder(structure_file, p1, p2):
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
    with open(script_path, 'r') as file:
        lines = file.readlines()
    
    with open(script_path, 'w') as file:
        for line in lines:
            if line.startswith("temp_dir="):
                file.write(f"temp_dir={dir_name}\n")
            else:
                file.write(line)


def main():
    """ Welcome """
    print("*" * 31)
    print()
    print((" " * 7) + "ColabFold Wrapper")
    print()
    print("*" * 31)

    jobs = "jobs.json"
    
    script_path = initialize_project(jobs)

    print(">>> ATTEMPTING TO RUN COLABFOLD")
    run_colabfold(script_path, jobs)

if __name__ == '__main__':
    main()