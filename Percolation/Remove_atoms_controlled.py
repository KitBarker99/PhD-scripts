from pymatgen.io.vasp.inputs import Poscar
import random

# Set a random seed for reproducibility
random_seed = 42
random.seed(random_seed)

# Load the input structure from POSCAR file
poscar_file = "all_1_layer_Zr.vasp"
structure = Poscar.from_file(poscar_file).structure

# Define the maximum distance for immediate neighbors and further radius for void check
max_distance = 3.0  # Radius for immediate neighbor check
further_radius = 6.0  # Radius to check for avoiding large voids

def remove_species_with_constraints(structure, chosen_species, removal_percentage, max_distance, further_radius):
    # Calculate the number of species to remove
    num_species = structure.composition[chosen_species]
    num_to_remove = int(num_species * (removal_percentage / 100))

    # Get indices of the chosen species
    species_indices = [i for i, site in enumerate(structure) if site.species_string == chosen_species]
    random.shuffle(species_indices)

    removed_atoms = 0

    while removed_atoms < num_to_remove and species_indices:
        for index in species_indices[:]:
            # Check distances to nearby atoms within max_distance
            site = structure[index]
            neighbors_within_max = structure.get_neighbors(site, max_distance)
            neighbors_within_further = structure.get_neighbors(site, further_radius)

            # Number of nearest neighbors within max_distance
            num_neighbors_within_max = len(neighbors_within_max)
            # Number of neighbors within further_radius
            num_neighbors_within_further = len(neighbors_within_further)

            # Check if there are neighbors within max_distance and that there are more than one neighbor within further_radius
            if any(neighbor.species_string != chosen_species for neighbor, _, _, _ in neighbors_within_max) and num_neighbors_within_further > 1:
                structure.remove_sites([index])
                removed_atoms += 1
                species_indices.remove(index)  # Remove index from the list after removal
                # Update species_indices for the new structure
                species_indices = [i for i, site in enumerate(structure) if site.species_string == chosen_species]
                random.shuffle(species_indices)  # Shuffle again after updating
                print(f"  --> Removed {chosen_species} atom at index {index}. Total removed: {removed_atoms}")
                print(f"      - Neighbors within {max_distance} Å: {num_neighbors_within_max}")
                print(f"      - Neighbors within {further_radius} Å: {num_neighbors_within_further}")
                break  # Reevaluate the structure after each removal
            else:
                print(f"  --> Did not remove {chosen_species} atom at index {index} (doesn't meet constraints)")
                print(f"      - Neighbors within {max_distance} Å: {num_neighbors_within_max}")
                print(f"      - Neighbors within {further_radius} Å: {num_neighbors_within_further}")

        # Exit loop if no more atoms can be removed
        if removed_atoms == num_to_remove or not species_indices:
            break

    print(f"\nSuccessfully removed {removed_atoms} {chosen_species} atoms.")
    return structure

# Remove species in a loop until the user decides to stop
while True:
    chosen_species = input("Choose a species to remove: ").strip()
    removal_percentage = float(input(f"Enter the removal percentage for {chosen_species}: "))

    structure = remove_species_with_constraints(structure, chosen_species, removal_percentage, max_distance, further_radius)
    
    more_removals = input("Do you want to remove another species? (yes/no): ").strip().lower()
    if more_removals != 'yes':
        break

# Save the modified structure to a new POSCAR file
new_poscar_file = f"final_structure_{structure.composition.reduced_formula}_random.vasp"
Poscar(structure).write_file(new_poscar_file)

# Print the reduced formula of the final structure
reduced_formula = structure.composition.reduced_formula
print("Species removal completed and saved to", new_poscar_file)
print("Reduced formula of the final structure:", reduced_formula)
print("Random seed used:", random_seed)