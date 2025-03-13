from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
import random

# Set a random seed for reproducibility
random_seed = 42
random.seed(random_seed)

# Load the input structure from POSCAR file
poscar_file = "Li2ZrCl6_5x5x10_all.vasp"
structure = Poscar.from_file(poscar_file).structure

# Define the maximum distance for immediate neighbors and further radius for void check
max_distance = 3.0  # Radius for immediate neighbor check
further_radius = 6.0  # Radius to check for avoiding large voids

def remove_species_simultaneously(structure):
    # Prompt user to choose species and specify the removal percentage for each
    chosen_species = input("Choose species to remove (comma-separated): ").split(',')
    removal_percentages = list(map(float, input("Enter removal percentages for each species (comma-separated): ").split(',')))

    if len(chosen_species) != len(removal_percentages):
        raise ValueError("The number of chosen species and removal percentages must match.")

    indices_to_remove = []

    for species, percentage in zip(chosen_species, removal_percentages):
        # Calculate the number of species to remove
        num_species = structure.composition[species]
        num_to_remove = int(num_species * (percentage / 100))

        # Get indices of the chosen species
        species_indices = [i for i, site in enumerate(structure) if site.species_string == species]
        random.shuffle(species_indices)

        removed_count = 0

        while removed_count < num_to_remove and species_indices:
            index = species_indices.pop()

            # Get immediate neighbors within max_distance
            immediate_neighbors = structure.get_neighbors(structure[index], max_distance)

            # Ensure removing this atom won't create a large void
            further_neighbors = structure.get_neighbors(structure[index], further_radius)

            # Print information about the removal process
            print(f"\nConsidering removal of {species} atom at index {index} with neighbors:")

            for neighbor in immediate_neighbors:
                neighbor_site = neighbor.site
                print(f" - Neighbor at index {structure.index(neighbor_site)} ({neighbor_site.species_string}), Distance: {neighbor.distance:.2f} Å")

            if len(further_neighbors) > 0:
                indices_to_remove.append(index)
                structure.remove_sites([index])  # Update structure immediately after removal
                print(f"Removal of {species} atom at index {index} accepted.\n")
                removed_count += 1
            else:
                print(f"Removal of {species} atom at index {index} rejected due to void check.\n")

    removed_species_counts = {species: indices_to_remove.count(i) for species in chosen_species}
    print(f"Removed species: {removed_species_counts}")

    return structure

# Remove species simultaneously based on user input
structure = remove_species_simultaneously(structure)

# Save the modified structure to a new POSCAR file
new_poscar_file = "Li2.8Zr0.8Cl6_modified.vasp"
Poscar(structure).write_file(new_poscar_file)

# Print the reduced formula of the final structure
reduced_formula = structure.composition.reduced_formula
print("Species removal completed and saved to", new_poscar_file)
print("Reduced formula of the final structure:", reduced_formula)
print("Random seed used:", random_seed)
