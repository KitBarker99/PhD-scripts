import os
import re
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.entries.mixing_scheme import MaterialsProjectDFTMixingScheme
from pymatgen.io.vasp import Vasprun
from pymatgen.entries.compatibility import ComputedEntry, ComputedStructureEntry
from pymatgen.core import Element

# Initialize API Key
API_KEY = ""

open_element = "Cu"
vasprun_file = "r2SCAN/vasprun.xml"
subdirs = sorted(set(next(os.walk('.'))[1]))  # Ensure no duplicates in subdirs
cwd = os.getcwd()  # get directory to return to

# Path to additional calculations
other_calculations_dir = os.path.abspath("../../../A2BMX6/copper/rubidium/Cl/")  # Ensure absolute path
other_calc_subdirs = next(os.walk(other_calculations_dir))[1]

def elements_to_sorted_string(elements_set):
    """
    Convert a set of Element objects into a sorted string representation.
    """
    element_symbols = sorted([el.symbol for el in elements_set])  # Extract symbols explicitly
    return "".join(element_symbols)

def extract_elements_from_formula(formula):
    """
    Extract element symbols from a chemical formula string and return as a sorted string.
    """
    return "".join(sorted(re.findall(r"[A-Z][a-z]?", formula)))

def find_matching_additional_entry(elements_set):
    """
    Find a suitable additional calculation that contains the same element set.
    """
    target_elements_str = elements_to_sorted_string(elements_set).strip()
    print(f"Searching for additional calculation with elements: {target_elements_str}")

    for folder in other_calc_subdirs:
        try:
            folder_elements_str = extract_elements_from_formula(folder).strip()
            print(f"Checking folder: {folder} with elements {folder_elements_str}")

            # Use strict string comparison after stripping
            if target_elements_str == folder_elements_str:
                print(f"Match found! Using {folder}")
                additional_vasprun_path = os.path.join(other_calculations_dir, folder, "r2SCAN", "vasprun.xml")
                additional_vasprun_path = os.path.abspath(additional_vasprun_path)  # Ensure absolute path
                print(f"Using vasprun with path: {additional_vasprun_path}")
                if os.path.exists(additional_vasprun_path):
                    print(f"Match confirmed! Using {folder} with path {additional_vasprun_path}")
                    return additional_vasprun_path  # Exit loop once match is found
                else:
                    print(f"File not found at {additional_vasprun_path}, continuing search.")
        except Exception as e:
            print(f"Error processing folder {folder}: {e}")

    print("Loop complete, no matching folder found.")
    return None

def create_adjusted_entry(entry, ehull, adjustment=0.1):
    """
    Creates a new ComputedStructureEntry or ComputedEntry with the adjusted energy.
    """
    original_energy = entry.energy
    adjusted_energy = original_energy - (ehull * entry.composition.num_atoms) - adjustment
    print(f"Original energy: {original_energy:.6f} eV")
    print(f"Adjusted energy: {adjusted_energy:.6f} eV")

    if isinstance(entry, ComputedStructureEntry):
        return ComputedStructureEntry(
            structure=entry.structure,
            energy=adjusted_energy,
            composition=entry.composition,
            energy_adjustments=entry.energy_adjustments,
            parameters=entry.parameters,
            data=entry.data,
            entry_id=entry.entry_id
        )
    else:
        return ComputedEntry(
            entry.composition,
            adjusted_energy,
            parameters=entry.parameters,
            data=entry.data,
            entry_id=entry.entry_id
        )

processed_subdirs = set()  # To track processed subdirectories

with open("new_r2SCAN.txt", "w") as sw_file, open("new_reactions.txt", "w") as reactions_file:
    sw_file.write("Structure\tLower limit (V)\tUpper limit (V)\tEnergy Above Hull (eV/atom)\tStability Window Error\tHull Calculation Error\n")

    for subdir in subdirs:
        if subdir in processed_subdirs:
            continue  # Skip if the subdir was already processed
        processed_subdirs.add(subdir)

        os.chdir(subdir)
        try:
            # Check if the vasprun.xml file exists
            if not os.path.exists(vasprun_file):
                print(f"Missing {vasprun_file} in {subdir}. Skipping...")
                sw_file.write(f"{subdir}\tfailed\tfailed\tfailed\n")
                os.chdir(cwd)
                continue

            vasprun = Vasprun(vasprun_file)
            entry = vasprun.get_computed_entry(inc_structure=True)
            if entry is None:
                sw_file.write(f"{subdir}\tfailed\tfailed\tfailed\n")
                os.chdir(cwd)
                continue

            # Check run_type
            if entry.parameters.get("run_type") == "r2SCAN":
                entry.parameters["run_type"] = "R2SCAN"
            elif entry.parameters.get("run_type") != "R2SCAN":
                print("No r2SCAN run found")  # Keep the print for debugging

            system_name = entry.composition.reduced_formula
            elements_set = set(entry.composition.elements)
            elements = set(vasprun.atomic_symbols)
            print(f"Processing {system_name} with elements: {elements_set}")

            additional_vasprun_path = find_matching_additional_entry(elements_set)
            additional_entry = None
            if additional_vasprun_path:
                additional_vasprun = Vasprun(additional_vasprun_path)
                additional_entry = additional_vasprun.get_computed_entry(inc_structure=True)
                if additional_entry.parameters.get("run_type") == "r2SCAN":
                    additional_entry.parameters["run_type"] = "R2SCAN"

            # Get unique entries from Materials Project
            with MPRester(API_KEY) as rester:
                mp_entries = rester.get_entries_in_chemsys(
                    elements=list(elements),
                    additional_criteria={"thermo_types": ["R2SCAN"]}
                )
                for mp_entry in mp_entries:
                    run_type = mp_entry.parameters.get("run_type")
                    if run_type == "r2SCAN":
                        mp_entry.parameters["run_type"] = "R2SCAN"
                        
                # Remove duplicates by entry_id
                mp_entries = list({entry.entry_id: entry for entry in mp_entries}.values())

            all_entries = [entry] + mp_entries
            if additional_entry and additional_entry not in all_entries:
                all_entries.append(additional_entry)
                print(f"Using additional entry for {system_name} from {additional_vasprun_path}")

            scheme = MaterialsProjectDFTMixingScheme()
            corrected_entries = scheme.process_entries(all_entries)
            pd = PhaseDiagram(corrected_entries)
            hull_calculation_error = "None"  # Default value in case no error occurs
            try:
                original_ehull = pd.get_e_above_hull(entry)
                print(f"The original energy above hull of {system_name} is {original_ehull:.3f} eV/atom.")
                if additional_entry is not None:
                    ehull_additional = pd.get_e_above_hull(additional_entry)
                    print(f"The original energy above hull of additional entry is {ehull_additional:.3f} eV/atom.")

                adjusted_entry = entry
                if original_ehull is not None and original_ehull > 0:
                    # Create an adjusted entry using the ComputedStructureEntry if needed
                    adjusted_entry = create_adjusted_entry(entry, original_ehull)
                    if additional_entry is not None:
                        corrected_entries = scheme.process_entries([adjusted_entry] + [additional_entry] + mp_entries)
                        pd = PhaseDiagram(corrected_entries)
                    else:
                        corrected_entries = scheme.process_entries([adjusted_entry] + mp_entries)
                        pd = PhaseDiagram(corrected_entries)
                        print(f"Using single entry without additional calculation for {subdir}")
                        sw_file.write(f"{system_name}\t{original_ehull:.3f}\tInfo: Only {subdir} entry used, additional entry was None.\n")
                else:
                    pd = PhaseDiagram(corrected_entries)
            except Exception as e:
                original_ehull = None
                hull_calculation_error = str(e)
                print(f"Failed to calculate energy above hull for {system_name}: {hull_calculation_error}")

            try:
                print(f"Current electrochemistry attempt using {subdir}")
                Li_entries = [e for e in corrected_entries if e.composition.reduced_formula == open_element]
                uLi0 = min(Li_entries, key=lambda e: e.energy_per_atom).energy_per_atom

                el_profile = pd.get_element_profile(Element(open_element), adjusted_entry.composition)

                voltages = []
                reactions_data = []
                for i, d in enumerate(el_profile):
                    voltage = -(d["chempot"] - uLi0)
                    voltages.append(voltage)
                    reactions_data.append((voltage, d["reaction"]))

                    stable_ranges = []
                    start_voltage = None
                    is_stable = lambda reactants, products: system_name in [r.reduced_formula for r in reactants] and \
                        system_name in [p.reduced_formula for p in products]

                for i, d in enumerate(el_profile):
                    if is_stable(d["reaction"].reactants, d["reaction"].products):
                        if start_voltage is None:
                            start_voltage = voltages[i]
                    else:
                        if start_voltage is not None:
                            stable_ranges.append((start_voltage, voltages[i]))
                            start_voltage = None
    
        # Ensure the last stable range is captured
                if start_voltage is not None:
                    stable_ranges.append((start_voltage, voltages[-1]))

                # Now, write results only ONCE
                if stable_ranges:
                    lower_limit = min(start for start, end in stable_ranges)
                    upper_limit = max(end for start, end in stable_ranges)
                    ehull_additional_str = f"{ehull_additional:.3f}" if additional_entry is not None else "None"
                    sw_file.write(f"{system_name}\t{lower_limit:.3f}\t{upper_limit:.3f}\t{original_ehull if original_ehull is not None else 'failed'}\t{ehull_additional_str}\t{hull_calculation_error}\n")
                else:
                    stability_window_error = "no_stable_ranges"
                    sw_file.write(f"{system_name}\tfailed\tfailed\t{original_ehull if original_ehull is not None else 'failed'}\t{stability_window_error}\t{hull_calculation_error}\n")

                    # Write reaction data only ONCE per system
                    reactions_file.write(f"Reactions for {system_name}:\n")
                for voltage, reaction in reactions_data:
                    reactions_file.write(f"Voltage: {voltage:.3f} V\n")
                    reactions_file.write(f"{reaction}\n\n")
                    
            except Exception as e:
                stability_window_error = str(e)
                sw_file.write(f"{system_name}\tfailed\tfailed\t{original_ehull if original_ehull is not None else 'failed'}\t{ehull_additional_str}\t{stability_window_error}\t{hull_calculation_error}\n")

        finally:
            os.chdir(cwd)
