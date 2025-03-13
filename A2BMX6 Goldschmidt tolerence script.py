import math

# Function to calculate the tolerance factor
def calculate_tolerance_factor(radii):
    r_A, r_B, r_M, r_X = radii
    t = (math.sqrt(2) * (r_A + r_X)) /  (r_B + r_M + (2 * r_X))
    return t

# List of elements and their ionic radii
elements = {
    'A': [('Rb', 1.72), ('Cs', 1.88), ('Tl', 1.7)],
    'B': [('Li', 0.76)],
    'M': [('Al', 0.535), ('Bi', 1.03), ('Ce', 1.01), ('Dy', 0.912), ('Er', 0.89),
          ('Eu', 0.947), ('Ga', 0.62), ('Gd', 0.938), ('Ho', 0.901), ('In', 0.80),
          ('La', 1.03), ('Lu', 0.861), ('Nd', 0.983), ('Pr', 0.99), ('Sc', 0.745),
          ('Sm', 0.958), ('Tm', 0.88), ('Y', 0.90), ('Yb', 0.868)],
    'X': [('Cl', 1.81), ('Br', 1.96), ('I', 2.2)]
}

# Open output file
with open("A2BMX6_tolerance_factors.txt", "w") as f_out:
    # Loop over all combinations of elements
    for (A_name, A_radii) in elements['A']:
        for (B_name, B_radii) in elements['B']:
            for (M_name, M_radii) in elements['M']:
                for (X_name, X_radii) in elements['X']:
                    radii = (A_radii, B_radii, M_radii, X_radii)
                    combination = f"{A_name}2{B_name}{M_name}{X_name}6"
                    t = calculate_tolerance_factor(radii)
                    # Write to output file
                    f_out.write(f"{combination} = {t:.4f}\n")
