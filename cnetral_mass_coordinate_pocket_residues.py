from Bio import PDB
import numpy as np
import os

def calculate_central_mass(residue):
    residue_atoms = [atom for atom in residue.get_atoms()]
    if len(residue_atoms) == 0:
        return None

    # Calculate the central mass of the residue
    residue_coordinates = [atom.get_coord() for atom in residue_atoms]
    central_mass = np.mean(residue_coordinates, axis=0)
    return central_mass

# Specify the path to your PDB file
refined_set_folder=os.listdir("refined-set")
for pdb_folder in refined_set_folder:
    pdb_file = "refined-set/"+pdb_folder+"/"+pdb_folder+"_pocket.pdb"

    # Create a PDB parser and parse the file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # Calculate and print the central mass coordinate for each residue
    for model in structure:
        for chain in model:
            for residue in chain:
                central_mass = calculate_central_mass(residue)
                if central_mass is not None:
                    print(f"Residue {residue.get_id()}: Central Mass Coordinate = {central_mass}")
~                                                                                                             
