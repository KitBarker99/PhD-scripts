from pymatgen.io.vasp.inputs import Poscar
import random

# Set a random seed for reproducibility
random_seed = 42
random.seed(random_seed)

# Load the input structure from POSCAR file
poscar_file = "all_1_layer_Zr.vasp"
structure = Poscar.from_file(poscar_file).structure

def remove_species(structure):
    # Prompt user to choose a species and specify the removal percentage
    chosen_species = input("Choose a species to remove: ")
    removal_percentage = float(input("Enter the removal percentage: "))

    # Calculate the number of species to remove
    num_species = structure.composition[chosen_species]
    num_to_remove = int(num_species * (removal_percentage / 100))

    # Get indices of the chosen species
    species_indices = [i for i, site in enumerate(structure) if site.species_string == chosen_species]
    random.shuffle(species_indices)

    # Remove the specified percentage of the chosen species
    indices_to_remove = species_indices[:num_to_remove]
    structure.remove_sites(indices_to_remove)

    print(f"Removed {num_to_remove} atoms of species {chosen_species}.")

    return structure

# Remove species in a loop until the user decides to stop
while True:
    structure = remove_species(structure)
    more_removals = input("Do you want to remove another species? (yes/no): ").strip().lower()
    if more_removals != 'yes':
        break

# Save the modified structure to a new POSCAR file
new_poscar_file = "0.66_Zr.vasp"
Poscar(structure).write_file(new_poscar_file)

# Print the reduced formula of the final structure
reduced_formula = structure.composition.reduced_formula
print("Species removal completed and saved to", new_poscar_file)
print("Reduced formula of the final structure:", reduced_formula)
print("Random seed used:", random_seed)
