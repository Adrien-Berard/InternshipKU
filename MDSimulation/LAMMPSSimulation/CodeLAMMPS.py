from lammps import lammps

def main():
    lmp = lammps()
    lmp.file("in.lj")

if __name__ == '__main__':
    main()