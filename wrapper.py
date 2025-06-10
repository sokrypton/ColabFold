"""
Luke Cirne
ColabFold Wrapper
Template Cycling with Experimental Distance Restraint Data
Ma Lab
"""
import os
import sys
import json
import getpass

def initialize_project(jobs):
    try:
        username = getpass.getuser()
        print(f">>> LOGGING AS: {username}")
    except OSError:
        print(">>> COULD NOT OBTAIN USERNAME")
        print(">>> RESULTS STORED UNDER DEFAULT")
        username = "DEFUALT"

    try:
        with open(jobs, "r") as file:
            print(">>> READING JSON")
            job_dict = json.load(file)
            user_JIDs = [int(job["JID"]) for job in job_dict["jobs"] if job["user"] == username]
            last_user_JID = max(user_JIDs) if user_JIDs else 0
            new_job = {
                "user": username,
                "JID": last_user_JID + 1
            }
            job_dict["jobs"].append(new_job)
        
        with open(jobs, 'w') as f:
            print(">>> APPENDING TO JSON")
            json.dump(job_dict, f, indent=4)
            print("APPEND COMPLETE")

    except FileNotFoundError:
        with open(jobs, "w") as file:
            print(">>> WRITING NEW JSON...")
            new_json = '''
                {
                    "jobs": [
                        {
                            "user": %(username)s,
                            "JID": 0
                        }
                    ]
                }
            ''' % {"username": username}
            print("WRITE COMPLETE")


def main():
    """ Welcome """
    print("*" * 31)
    print()
    print((" " * 7) + "ColabFold Wrapper")
    print()
    print("*" * 31)




if __name__ == '__main__':
    main()