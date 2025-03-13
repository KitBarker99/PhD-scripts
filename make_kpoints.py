from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

struct = Structure.from_file("test.cif")
prim_struct=SpacegroupAnalyzer(struct)
new_struct=prim_struct.get_primitive_standard_structure(international_monoclinic=True)
kpath = HighSymmKpath(new_struct)
kpts = Kpoints.automatic_linemode(divisions=40,ibz=kpath)
kpts.write_file("KPOINTS_nsc")
