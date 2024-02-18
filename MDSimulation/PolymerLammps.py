from lammps import PyLammps
L = PyLammps()
L.command("region box block 0 10 0 5 -0.5 0.5")
L.system()
L.system.xlo, L.system.xhi = 0.0, 20.0
L.system.ylo, L.system.yhi = 0.0, 20.0
L.system.zlo, L.system.zhi = 0.0, 20.0