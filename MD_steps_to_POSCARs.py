#   This script is designed to read VASP output (XDATCAR) from a MD calculation
#   and subsequently produce a number of POSCAR files at different steps
#   for general minima search via relaxation of the positionsat these steps.



#READ BEFORE USE

#   Depending on size and length of the MD run this may result in a large
#   number of files and folders i.e. 200 are created in this example

#   Please analyse XDATCAR before use and choose 'split' appropriately.
    

#   Import required modules

import os
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet

#------------------------------------------------------------------------

path = os.getcwd()

# read the lines
f=open("XDATCAR")
lines=f.readlines()
f.close()


no_atoms = input("Number of atoms in the structure: ")
split = input("Number of steps between each POSCAR generation: ")
atoms = int(no_atoms)
steps = int(split)

file_number=0 # counter for the files
line_counter=0 # keeps track of the number of lines written in the file so far

for i in range(7,len(lines)):
    if line_counter==0: # we started writing to a file
        os.mkdir(path + "/step" + str(file_number)) # make new directory
        os.chdir(path + "/step" + str(file_number)) # move to new directory
        f=open("POSCAR","w") # open file
        for j in range(7):# write the first seven lines
            f.write(lines[j])

    if line_counter <= atoms:
        f.write(lines[i])# write one line
        line_counter+=1 # increment line counter
    else:
        line_counter+=1
    
    if line_counter == (atoms + 1) * steps:
        f.close()
        
#------------------------------------------------------------------------------------------
#       Following code creates input files for a VASP relaxation using MP parameters
        
        Structure = Structure.from_file("POSCAR")
        relax = MPRelaxSet(Structure,user_incar_settings={"ISMEAR":0,"EDIFF":1E-5,"ISYM":0,"NCORE":4,"NSW":191})
        relax.write_input("structure_relax")
        
        os.chdir("structure_relax")

#-------------------------------------------------------------------------------
#   Creates a script for the HPC to run VASP, change bottom line appropriately

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
        
#-------------------------------------------------------------------------------------

        
        file_number+=1 # start writing to a new file
        line_counter=0 
        os.chdir(path)

if f.closed == False: # the file was not closed properly in the for loop
    f.close()
    
#---------------------------------------------------------------------
#       Done  