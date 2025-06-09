"""
Luke Cirne
ColabFold Wrapper
Template Cycling with Experimental Distance Restraint Data
Ma Lab
"""
import os
import sys
import json

def initialize_project(jobs):
    try:
        with open(jobs, "r") as file:
            print("reading")
            job_dict = json.load(file)
            
    except FileNotFoundError:
        with open(jobs, "w") as file:
            print("writing")


def main():
    """ Welcome """
    print("*" * 31)
    print()
    print((" " * 7) + "ColabFold Wrapper")
    print()
    print("*" * 31)


if __name__ == '__main__':
    main()