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

def compute_distances(atoms):
    pairs_list = ["AA", "AU", "AC", "AG", "UU", "UC", "UG", "CC", "CG", "GG"]
    pairs = {pair: {"distances": []} for pair in pairs_list}
    pairs["all_distances"] = []

    for key in atoms:
        for key2 in atoms:
            if key2 > key:
                if atoms[key2]["position"] - atoms[key]["position"] > 3:
                    dist = int(np.sqrt((atoms[key2]['x'] - atoms[key]['x'])**2 +
                                       (atoms[key2]['y'] - atoms[key]['y'])**2 +
                                       (atoms[key2]['z'] - atoms[key]['z'])**2))

                    if dist <= 20:
                        pair = atoms[key]["nucleotide"] + atoms[key2]["nucleotide"]
                        if pair in pairs:
                            pairs[pair]['distances'].append(dist)
                        elif pair[::-1] in pairs:  # Check the reverse order of nucleotides
                            pairs[pair[::-1]]['distances'].append(dist)
                        pairs["all_distances"].append(dist)

    return pairs

def compute_frequencies(pairs):
    pairs["frequencies_ref"] = [pairs['all_distances'].count(i) / len(pairs['all_distances']) for i in range(20)]

    for pair in pairs:
        if pair != "all_distances" and pair != "frequencies_ref":
            pairs[pair]['frequencies'] = [pairs[pair]['distances'].count(i) / len(pairs[pair]['distances']) for i in range(20)]

def compute_scores(pairs):
    for pair in pairs:
        if pair != "all_distances" and pair != "frequencies_ref":
            pairs[pair]["scores"] = [min(10, -log(pairs[pair]['frequencies'][i] / pairs['frequencies_ref'][i]))
                                     if pairs['frequencies_ref'][i] != 0 and pairs[pair]['frequencies'][i] != 0
                                     else 0 for i in range(20)]

def save_scores_to_files(pairs):
    for pair in pairs:
        if pair != "all_distances" and pair != "frequencies_ref":
            filename = f"data/{pair}.txt"
            with open(filename, 'w') as f:
                f.writelines([f"{score}\n" for score in pairs[pair]['scores']])

def main():
    # Get the input file as an argument
    input_file = "DATA/8bu8.pdb"

    # Read PDB file and create atoms dictionary
    pdb_data = read_pdb_file(input_file)
    atoms = create_atoms_dict(pdb_data)

    # Compute distances between residues
    pairs = compute_distances(atoms)

    # Compute frequencies
    compute_frequencies(pairs)

    # Compute scores
    compute_scores(pairs)

    # Save scoring data into files
    save_scores_to_files(pairs)

if __name__ == "__main__":
    main()
