from lammps import lammps, PyLammps
import numpy as np

natoms = 1000
system_size = 20.0 # In Angstrom

positions = []
for i in range(natoms):
    positions.append(np.random.rand(3)*system_size)

# Write LAMMPS data file
with open("random.data",'w') as fdata:
    fdata.write('Random atoms - written for EnvodeVentorTutorial\n\n')

    # --- HEADER ---#
    # Specify number of atoms and atom types
    fdata.write('{} atoms\n'.format(natoms))
    fdata.write('{} atom types\n'.format(3))
    # specify  box dimesions
    fdata.write('{} {} xlo xhi\n'.format(0.0, system_size))
    fdata.write('{} {} ylo yhi\n'.format(0.0, system_size))
    fdata.write('{} {} zlo zhi\n'.format(0.0, system_size))

    # Atoms section
    fdata.write('Atoms\n\n')

    # Write each position
    for pos in enumerate(positions):
        state = np.random.choice([1,2,3])
        fdata.write('{} '+str(state)+' {} {} {}\n'.format(i+1,*pos))

    # Write bonds