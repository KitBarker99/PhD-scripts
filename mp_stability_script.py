#!/usr/local/bin/python
"""
This is a basic example of how to create, plot, and analyze OPEN Phase Diagrams using the pymatgen
codebase and Materials Project database. To run this example, you should:
* have pymatgen (www.pymatgen.org) installed along with matplotlib
* obtain a Materials Project API key (https://www.materialsproject.org/open)
* paste that API key in the MAPI_KEY variable below, e.g. MAPI_KEY = "foobar1234"
For citation, see https://www.materialsproject.org/citing
For the accompanying comic book, see http://www.hackingmaterials.com/pdcomic
"""
from __future__ import print_function

import pymatgen
from pymatgen import io
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.ext.matproj import MPRester, Element,Composition
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram, \
    PhaseDiagram, PDPlotter

if __name__ == "__main__":

    mpr = MPRester("") 
    vasprun = pymatgen.io.vasp.outputs.Vasprun("C:/Users/kitba/OneDrive/Documents/PhD/CIFs/Li3MX6/La/trig/vasprun.xml")
    my_entry = vasprun.get_computed_entry(inc_structure=False)
    print(mpr.get_stability([my_entry])[0])
