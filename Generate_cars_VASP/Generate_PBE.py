import os
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet
import shutil
from pymatgen.io.vasp import Kpoints
from pymatgen.io.vasp.sets import MPScanRelaxSet


root_directory = os.getcwd()

#uses a list of elements to substitute into a POSCAR
with open(elements_file, "r") as f:
    elements_list = [line.strip() for line in f]

# Read the original POSCAR structure
original_structure = Poscar.from_file(poscar_file).structure

# Loop over each replacement element
for new_element in elements_list:
    # Create a copy of the original structure
    modified_structure = original_structure.copy()

    # Replace one of the elements with the new elements
    for i, site in enumerate(modified_structure):
        if site.species_string == "In":
            modified_structure.replace(i, new_element)

    # Create a new directory for each substitution
    # Name accordingly
    output_directory = f"Cs3Li4{new_element}2Cl13"
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Write the modified POSCAR file to the new directory
    poscar_filename = os.path.join(output_directory, "POSCAR")
    Poscar(modified_structure).write_file(poscar_filename)
    os.chdir(root_directory)

print("Substitution complete. New POSCAR files saved in the appropriate directories.")

# Ask the user for the type of calculation
calculation_type = input("Enter the type of calculation to generate (1 for PBE, 2 for r2SCAN): ")

if calculation_type == '1':
    # PBE calculation
    for directory_name in os.listdir():
        if os.path.isdir(directory_name) and 'structure_relax' not in os.listdir(directory_name):
            if 'POSCAR' in os.listdir(directory_name):
                # Change to the current subdirectory
                os.chdir(directory_name)
                # Run the script
                structure = Structure.from_file("POSCAR")
                relax = MPRelaxSet(structure, user_incar_settings={"NELM":300,"ISMEAR":0,"EDIFF":1E-5,"ISYM":0,"NCORE":4,"NSW":191})
                relax.write_input("structure_relax")
                os.chdir(root_directory)

elif calculation_type == '2':
    # r2SCAN calculation
    # Set the k-point grid size using Monkhorst-Pack
    kpoints = Kpoints.monkhorst_automatic(kpts=(2, 2, 2), shift=(0, 0, 0))

    for directory_name in os.listdir():
        if os.path.isdir(directory_name) and 'r2SCAN_relax' not in os.listdir(directory_name):
            if 'POSCAR' in os.listdir(directory_name):
                # Change to the current subdirectory
                os.chdir(directory_name)
                # Run the script
                structure = Structure.from_file("POSCAR")
                relax = MPScanRelaxSet(structure, auto_kspacing=False, user_kpoints_settings=kpoints, user_incar_settings={"ENCUT":600,"NELM":200,"LCHARG":"TRUE","ISMEAR":0,"LWAVE":"TRUE","EDIFF":1E-5,"ISYM":0,"NCORE":4,"NSW":191})
                relax.write_input("r2SCAN_relax")
                os.chdir(root_directory)

else:
    print("Error. Please input one of the specified values.")