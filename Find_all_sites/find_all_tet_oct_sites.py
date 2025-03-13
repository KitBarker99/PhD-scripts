from itertools import combinations
from pymatgen.core import Structure, Element
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.cif import CifWriter
import numpy as np

structure = Structure.from_file("Cu_missing.vasp")
cutoff = float(input("Type the cutoff Radius:")) # Define the cutoff distance


# Find the species in the structure
species = [site.species.reduced_formula for site in structure]

# Determine the target anion from the species
anion_list = ["F2", "Cl2", "Br", "I2", "O2", "S"]
anion = None
for sp in set(species):
    if sp in anion_list:
        anion = sp
        break

if anion is None:
    raise ValueError("No anion found in the structure")

anion_sites = []

# Iterate over each site in the structure
for site in structure:
    if site.species.reduced_formula == anion:
        anion_sites.append(site)

#print(anion_sites)



tetrahedra_centers = []
octahedra_centers = []

# Iterate over each anion site
for site in anion_sites:
    anions_in_polyhedra = []
    for candidate_site in anion_sites:
        if candidate_site == site:
            continue
        distance = np.linalg.norm(site.coords - candidate_site.coords)
        if distance < cutoff:
            anions_in_polyhedra.append(candidate_site)
    for n_anions in range(3, 6):
        for combination in combinations(anions_in_polyhedra, n_anions):
            if n_anions == 3:
                center = (site.coords + np.sum([a.coords for a in combination], axis=0)) / (n_anions + 1)
                tetrahedra = [site.coords, *[a.coords for a in combination]]
                angles = []
                for i in range(4):
                    for j in range(i + 1, 4):
                        for k in range(j + 1, 4):
                            i_vec = tetrahedra[i] - center
                            j_vec = tetrahedra[j] - center
                            k_vec = tetrahedra[k] - center
                            angles.append(np.arccos(np.dot(i_vec, j_vec) / (np.linalg.norm(i_vec) * np.linalg.norm(j_vec))))
                            angles.append(np.arccos(np.dot(j_vec, k_vec) / (np.linalg.norm(j_vec) * np.linalg.norm(k_vec))))
                            angles.append(np.arccos(np.dot(k_vec, i_vec) / (np.linalg.norm(k_vec) * np.linalg.norm(i_vec))))
                #print("Tetrahedral angles:", np.degrees(angles))
                if all(106 < a < 113 for a in np.degrees(angles)):
                    tetrahedra_centers.append(center)
            elif n_anions == 5:
                center = (site.coords + np.sum([a.coords for a in combination], axis=0)) / (n_anions + 1)
                octahedra = [site.coords, *[a.coords for a in combination]]
                angles = []
                for i in range(4):
                    for j in range(i + 1, 4):
                        i_vec = octahedra[i] - center
                        j_vec = octahedra[j] - center
                        angles.append(np.arccos(np.dot(i_vec, j_vec) / (np.linalg.norm(i_vec) * np.linalg.norm(j_vec))))
                #print("Octahedral angles:", np.degrees(angles))
                if all((86.5 < a < 93.5) or (177 < a < 183) for a in np.degrees(angles)):
                    octahedra_centers.append(center)

#print("Tetrahedra centers:", tetrahedra_centers)
#print("Octahedra centers:", octahedra_centers)

tet_atom_symbol = input("Enter the symbol of the atom you want to add to the tetrahedral sites: ")
tet_element = Element(tet_atom_symbol)
tet_structure = Structure(lattice=structure.lattice, species=[tet_element] * len(tetrahedra_centers),
                        coords=tetrahedra_centers, coords_are_cartesian=True)

#tet_structure.to(filename="tet_structure.vasp", fmt="poscar")

oct_atom_symbol = input("Enter the symbol of the atom you want to add to the octahedral sites: ")
oct_element = Element(oct_atom_symbol)
oct_structure = Structure(lattice=structure.lattice, species=[oct_element] * len(octahedra_centers),
                        coords=octahedra_centers, coords_are_cartesian=True)

# Concatenate the coordinates of the original structure and the two new structures
new_sites = structure.sites + tet_structure.sites + oct_structure.sites
new_structure = Structure.from_sites(new_sites)

new_structure.to(filename="structure_with_added_atoms.vasp", fmt="poscar")

cif_writer = CifWriter(new_structure, symprec=0.1, angle_tolerance=10)
cif_writer.write_file("new_structure.cif")

