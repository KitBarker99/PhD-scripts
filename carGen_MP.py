
# NOTES:

# MUST DO 'SOURCE ACTIVATE my_pymatgen'
# Change user_incar_settings to desired parameters
# POSCAR required in directory where this script is being called from



from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet
Structure = Structure.from_file("POSCAR")
relax = MPRelaxSet(Structure,user_incar_settings={"ISMEAR":0,"EDIFF":1E-5,"ISYM":0,"NCORE":4,"NSW":191})
relax.write_input("structure_relax")





#===================================== VASP_Run_Script generator ====================================

import os

#path = os.getcwd()
os.chdir("structure_relax")

with open('vasp_script','w') as f:
    f.writelines('''
#PBS -lselect=1:ncpus=32:mpiprocs=32:mem=60gb
#PBS -lwalltime=48:00:00
#PBS -N vasprun

# Load modules for applications

module load mpi/intel-2019.8.254
module load intel-suite/2017.6
module load gcc

# Change to the directory the job was submitted from

cd $PBS_O_WORKDIR

# Run program, using 'mpiexec' to start the job
# mpiexec automatically picks up the # of cores
# assigned to the job. No other flags are required
#  - note: don't use 'mpirun'
# 'mpiexec /absolute path to vasp executable'

mpiexec /rds/general/user/kab121/home/VASP/vasp.6.1.2_patched_vtst/bin/vasp_std
''')
