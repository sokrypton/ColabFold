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

def initialize_project(jobs) -> list:
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
        key_values.append(("JID", 0))
    append_json(jobs, key_values)

    # Obtain input file
    while True:
        input_file = input("Input file name in the format 'name.fasta': ")
        if os.path.isfile(f"./{input_file}"):
            break
        else:
            print("###### Invalid filepath ######")
    
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
    

    # Obtain num seeds


    # Create shell script to run colabfold_batch
    current_dir = os.path.dirname(os.path.abspath(__file__))
    script_name = "wrapper.sh"
    script_path = os.path.join(current_dir, script_name)
    # Writing script
    script_content = f"""
        JID={current_JID}
        num_c=
    """
     

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


def main():
    """ Welcome """
    print("*" * 31)
    print()
    print((" " * 7) + "ColabFold Wrapper")
    print()
    print("*" * 31)




if __name__ == '__main__':
    main()