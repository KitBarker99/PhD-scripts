# Convert the list of cartesian coordinates to fractional coordinates using the lattice parameters in the structure
fractional_tet_coords = [np.linalg.solve(structure.lattice.matrix, cart_coord) for cart_coord in tetrahedra_centers]
fractional_oct_coords = [np.linalg.solve(structure.lattice.matrix, cart_coord) for cart_coord in octahedra_centers]


#print("Fractional tetrahedral coords:", fractional_tet_coords)
#print("Fractional octahedral coords:", fractional_oct_coords)


tet_site_species = str(input("\nType the species you would like on the tetrahedral sites: "))
oct_site_species = str(input("\nType the species you would like on the octahedral sites: "))


# Create a list of PeriodicSite objects with the new atomic species and coordinates
new_tet_sites = [PeriodicSite(tet_site_species, coord, structure.lattice) for coord in fractional_tet_coords]
new_oct_sites = [PeriodicSite(oct_site_species, coord, structure.lattice) for coord in fractional_oct_coords]


# Combine the new sites with the existing sites in the Structure object
all_sites = structure.sites + new_tet_sites + new_oct_sites

# Create a new Structure object with the modified list of sites
structure = Structure(structure.lattice, [site.species for site in all_sites], [site.frac_coords for site in all_sites])



# Append the new sites to the existing structure
#structure.extend(new_tet_sites)
#structure.extend(new_oct_sites)

# Save the new Structure object to a file in VASP format
structure.to(fmt="POSCAR", filename="all_positions.vasp")


print("\nCongratulations, your new structure has been created as 'all_positions.vasp'")

#cif_writer = CifWriter(structure)
#cif_writer.write_file("new_structure.cif")
