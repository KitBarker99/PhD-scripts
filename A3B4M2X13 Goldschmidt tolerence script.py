import math

# Function to calculate the tolerance factor
def calculate_tolerance_factor(r_A_eff, r_M, r_X):
    t = (r_A_eff + r_X) / (math.sqrt(2) * (r_M + r_X))
    return t

# List of elements and their ionic radii
elements = {
    'A': [('Rb', 1.72), ('Cs', 1.88), ('Tl', 1.7)],
    'B': [('Li', 2.34), ('Cu', 2.29)],
    'M': [('Al', 0.535), ('Bi', 1.03), ('Ce', 1.01), ('Dy', 0.912), ('Er', 0.89),
          ('Eu', 0.947), ('Ga', 0.62), ('Gd', 0.938), ('Ho', 0.901), ('In', 0.80),
          ('La', 1.03), ('Lu', 0.861), ('Nd', 0.983), ('Pr', 0.99), ('Sc', 0.745),
          ('Sm', 0.958), ('Tm', 0.88), ('Y', 0.90), ('Yb', 0.868)],
    'X': [('Cl', 1.81), ('Br', 1.96), ('I', 2.2)]
}

# Open output file
with open("A3B4M2X13_tolerance_factors.txt", "w") as f_out:
    # Loop over all combinations of elements
    for (A_name, A_radii) in elements['A']:
        for (B_name, B_radii) in elements['B']:
            for (M_name, M_radii) in elements['M']:
                for (X_name, X_radii) in elements['X']:
                    # Calculate r_A_eff using the given formula
                    r_A_eff = ((3/4) * A_radii) + ((1/4) * B_radii)
                    radii = (r_A_eff, M_radii, X_radii)
                    combination = f"{A_name}3{B_name}4{M_name}2{X_name}13"
                    t = calculate_tolerance_factor(*radii)
                    # Write to output file
                    f_out.write(f"{combination} = {t:.4f}\n")
