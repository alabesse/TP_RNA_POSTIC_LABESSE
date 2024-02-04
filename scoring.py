import os
import re
import numpy as np
from math import log
import sys

def read_pdb_file(input_file):
    with open(input_file, "r") as f:
        lines = f.readlines()
    return [line for line in lines if re.search("ATOM[\d\s\D]*C3'", line)]

def create_atoms_dict(data):
    atoms = {}
    atom_id = 0
    for line in data:
        fields = line.split()
        if fields[4] == 'A':
            atoms[atom_id] = {
                "nucleotide": fields[3],
                "chain_id": fields[4],
                "position": int(fields[5]),
                "x": float(fields[6]),
                "y": float(fields[7]),
                "z": float(fields[8])
            }
            atom_id += 1
    return atoms

def compute_distances_and_scores(atoms):
    distances = []

    for key in atoms:
        for key2 in atoms:
            if key2 > key and atoms[key2]["position"] - atoms[key]["position"] == 4:
                dist = int(np.sqrt((atoms[key2]['x'] - atoms[key]['x'])**2 +
                                   (atoms[key2]['y'] - atoms[key]['y'])**2 +
                                   (atoms[key2]['z'] - atoms[key]['z'])**2))

                if dist <= 20:
                    distances.append(dist)

    return distances

def linear_interpolation(x, x1, y1, x2, y2):
    return y1 + ((y2 - y1) / (x2 - x1)) * (x - x1)

def calculate_gibbs_free_energy(distances, scores):
    total_score = 0

    for dist in distances:
        if 0 <= dist < len(scores):
            # Linear interpolation of scores
            idx = min(dist, 20)
            score = linear_interpolation(dist, idx-1, scores[idx-1], idx, scores[idx])
            total_score += score

    return total_score

def main():
    # Get the input file as an argument
    input_file = "DATA/8bu8.pdb"

    # Read PDB file and create atoms dictionary
    pdb_data = read_pdb_file(input_file)
    atoms = create_atoms_dict(pdb_data)

    # Specify the folder where the scores files are stored
    scores_folder = "DATA"

    # Loop through all score files
    for score_file in os.listdir(scores_folder):
        if score_file.endswith(".txt"):
            scores = [float(line.strip()) for line in open(os.path.join(scores_folder, score_file), 'r')]

            # Compute distances and calculate estimated Gibbs free energy
            distances = compute_distances_and_scores(atoms)
            gibbs_free_energy = calculate_gibbs_free_energy(distances, scores)

            print(f"File: {score_file}, Estimated Gibbs Free Energy: {gibbs_free_energy}")

if __name__ == "__main__":
    main()
