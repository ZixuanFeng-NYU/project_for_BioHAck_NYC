import numpy as np
from Bio import PDB
import os

def parse_mol2(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atom_data = []
    reading_atoms = False

    for line in lines:
        if line.strip() == "@<TRIPOS>ATOM":
            reading_atoms = True
        elif line.strip() == "@<TRIPOS>BOND":
            reading_atoms = False
        elif reading_atoms:
            fields = line.split()
            atom_data.append({
                'x': float(fields[2]),
                'y': float(fields[3]),
                'z': float(fields[4])
            })

    return atom_data

def calculate_com(coordinates):
    if len(coordinates) == 0:
        return None

    com = np.mean(coordinates, axis=0)
    return com

def calculate_distance(point1, point2):
    return np.linalg.norm(point1 - point2)

# Specify the path to your PDB and MOL2 files
refined_set_folder = os.listdir("refined-set")
for pdb_folder in refined_set_folder:
    pocket_file = "refined-set/" + pdb_folder + "/" + pdb_folder + "_pocket.pdb"
    ligand_file = "refined-set/" + pdb_folder + "/" + pdb_folder + "_ligand.mol2"

    # Parse the ligand MOL2 file
    ligand_data = parse_mol2(ligand_file)
    ligand_com = calculate_com(np.array([[atom['x'], atom['y'], atom['z']] for atom in ligand_data]))

    # Create a PDB parser and parse the pocket PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pocket_file)

    # Create a dictionary to store residue central mass distances
    residue_distances = {}

    # Calculate and store the central mass coordinate for each residue
    for model in structure:
        for chain in model:
            for residue in chain:
                central_mass = calculate_com(np.array([atom.get_coord() for atom in residue.get_atoms()]))
                if central_mass is not None:
                    residue_distances[residue.get_id()] = calculate_distance(central_mass, ligand_com)

    # Sort the residues by distance to the ligand central mass
    sorted_residues = sorted(residue_distances.items(), key=lambda x: x[1])

    # Print the ranked residues
    print(f"Ranking of residues based on distance to ligand central mass in {pdb_folder}:")
    for rank, (residue_id, distance) in enumerate(sorted_residues, start=1):
        print(f"Rank {rank}: Residue {residue_id}, Distance = {distance}")
