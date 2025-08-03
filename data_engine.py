"""
Luke Cirne
Data Engine for ColabFold Wrapper
Functions for graphing and data analysis
Ma Lab
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def graph_output_accuracy(distances: dict) -> str:
    # Collect and convert distances
    dist_values = [float(d) for d in distances.keys()]

    # Define bin range: from 5 below min to 5 above max, in steps of 10
    min_d = min(dist_values)
    max_d = max(dist_values)
    bin_start = np.floor(min_d - 5)
    bin_end = np.ceil(max_d + 5)
    bin_edges = np.arange(bin_start, bin_end + 10, 10)  # +10 to include final edge

    # Plot
    plt.figure(figsize=(8, 5))
    plt.hist(dist_values, bins=bin_edges, edgecolor="black", color="skyblue")
    plt.title("CF Output Distances (Å)")
    plt.xlabel("Distance (Å)")
    plt.ylabel("Frequency")
    plt.xticks(bin_edges)
    plt.tight_layout()

    # Save
    plot_name = "iteration_distances_hist"
    plt.savefig(f"{plot_name}.png")
    return plot_name


def test_graph_output_accuracy():
    test_dict = {
        110: "file0",
        89: "file1",
        104: "file2"
    }
    graph_output_accuracy(test_dict)


def testing():
    test_graph_output_accuracy()

if __name__ == "__main__":
    testing()