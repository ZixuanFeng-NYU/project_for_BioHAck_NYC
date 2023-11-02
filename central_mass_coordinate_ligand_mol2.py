import numpy as np

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

def calculate_com_ligand(mol2_data):
    if len(mol2_data) == 0:
        return None

    coordinates = np.array([[atom['x'], atom['y'], atom['z']] for atom in mol2_data])
    com_ligand = np.mean(coordinates, axis=0)
    return com_ligand

# Specify the path to your MOL2 file
mol2_file = "refined-set/1ebw/1ebw_ligand.mol2"

# Parse the MOL2 file
mol2_data = parse_mol2(mol2_file)

com_ligand = calculate_com_ligand(mol2_data)

if com_ligand is not None:
    print(f"Center of Mass (COM) of Ligand: {com_ligand}")
else:
    print("No atom data found in the MOL2 file.")

