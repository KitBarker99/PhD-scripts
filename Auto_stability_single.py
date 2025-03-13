import os
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.entries.mixing_scheme import MaterialsProjectDFTMixingScheme
from pymatgen.io.vasp import Vasprun
from pymatgen.entries.compatibility import ComputedEntry, ComputedStructureEntry
from pymatgen.core import Element

# Initialize API Key
API_KEY = ""

open_element = "Li"
vasprun_file = "r2SCAN/vasprun.xml"
subdirs = next(os.walk('.'))[1]
subdirs.sort()  # Sort the subdirectories alphabetically

def create_adjusted_entry(entry, ehull, adjustment=0.1):
    """
    Creates a new ComputedStructureEntry or ComputedEntry with the adjusted energy.
    If the entry has a structure (ComputedStructureEntry), create a new adjusted version of that.
    """
    original_energy = entry.energy
    adjusted_energy = original_energy - (ehull * entry.composition.num_atoms) - adjustment

    # Print the original and adjusted energy
    print(f"Original energy: {original_energy:.6f} eV")
    print(f"Adjusted energy: {adjusted_energy:.6f} eV")

    # Create a new entry with adjusted energy
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

with open("stability_windows_r2SCAN.txt", "w") as sw_file, open("reactions.txt", "w") as reactions_file:
    sw_file.write("Structure\tLower limit (V)\tUpper limit (V)\tEnergy Above Hull (eV/atom)\tStability Window Error\tHull Calculation Error\n")

    for subdir in subdirs:
        os.chdir(subdir)
        stability_window_error = ""
        hull_calculation_error = ""
        try:
            vasprun = Vasprun(vasprun_file)
            entry = vasprun.get_computed_entry(inc_structure=True)
            if entry is None:
                stability_window_error = "entry_is_none"
                hull_calculation_error = "entry_is_none"
                sw_file.write(f"{subdir}\tfailed\tfailed\tfailed\t{stability_window_error}\t{hull_calculation_error}\n")
                os.chdir('..')
                continue

            system_name = entry.composition.reduced_formula

            # Check whether the run has converged (ionic and electronic)
            is_ionic_converged = vasprun.converged_ionic
            max_electronic_steps = vasprun.parameters["NELM"]
            is_electronic_converged = all(len(step["electronic_steps"]) < max_electronic_steps for step in vasprun.ionic_steps)

            if not is_ionic_converged or not is_electronic_converged:
                stability_window_error = "ionic_or_electronic_not_converged"
                hull_calculation_error = "ionic_or_electronic_not_converged"
                sw_file.write(f"{system_name}\tfailed\tfailed\tfailed\t{stability_window_error}\t{hull_calculation_error}\n")
            else:
                elements = set(vasprun.atomic_symbols)

                # Fetch entries using the new Materials Project API
                with MPRester(API_KEY) as rester:
                    mp_entries = rester.get_entries_in_chemsys(
                        elements=list(elements), 
                        additional_criteria={"thermo_types": ["R2SCAN"]}
                    )

                if not mp_entries:
                    sw_file.write(f"{system_name}\tfailed\tfailed\tfailed\tno_entries_found\t\n")
                    os.chdir('..')
                    continue

                # Apply corrections with the mixing scheme
                scheme = MaterialsProjectDFTMixingScheme()
                corrected_entries = scheme.process_entries([entry] + mp_entries)

                # Create phase diagram
                pd = PhaseDiagram(corrected_entries)
                plotter = PDPlotter(pd)

                try:
                    original_ehull = pd.get_e_above_hull(entry)
                    print(f"The original energy above hull of {system_name} is {original_ehull:.3f} eV/atom.")

                    adjusted_entry = entry
                    if original_ehull is not None and original_ehull > 0:
                        # Create an adjusted entry using the ComputedStructureEntry if needed
                        adjusted_entry = create_adjusted_entry(entry, original_ehull)
                        corrected_entries = scheme.process_entries([adjusted_entry] + mp_entries)
                        pd = PhaseDiagram(corrected_entries)
                    else:
                        pd = PhaseDiagram(corrected_entries)
                except Exception as e:
                    original_ehull = None
                    hull_calculation_error = str(e)
                    print(f"Failed to calculate energy above hull for {system_name}: {hull_calculation_error}")

                try:
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

                    if start_voltage is not None:
                        stable_ranges.append((start_voltage, voltages[-1]))

                    if stable_ranges:
                        lower_limit = min(start for start, end in stable_ranges)
                        upper_limit = max(end for start, end in stable_ranges)
                        sw_file.write(f"{system_name}\t{lower_limit:.3f}\t{upper_limit:.3f}\t{original_ehull if original_ehull is not None else 'failed'}\t\t{hull_calculation_error}\n")
                    else:
                        stability_window_error = "no_stable_ranges"
                        sw_file.write(f"{system_name}\tfailed\tfailed\t{original_ehull if original_ehull is not None else 'failed'}\t{stability_window_error}\t{hull_calculation_error}\n")

                    reactions_file.write(f"Reactions for {system_name}:\n")
                    for voltage, reaction in reactions_data:
                        reactions_file.write(f"Voltage: {voltage:.3f} V\n")
                        reactions_file.write(f"{reaction}\n\n")
                except Exception as e:
                    stability_window_error = str(e)
                    sw_file.write(f"{system_name}\tfailed\tfailed\t{original_ehull if original_ehull is not None else 'failed'}\t{stability_window_error}\t{hull_calculation_error}\n")
        except Exception as e:
            stability_window_error = str(e)
            hull_calculation_error = str(e)
            sw_file.write(f"{subdir}\tfailed\tfailed\tfailed\t{stability_window_error}\t{hull_calculation_error}\n")
        os.chdir('..')

