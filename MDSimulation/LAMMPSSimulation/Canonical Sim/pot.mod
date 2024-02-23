# Exponential potential
pair_style buck 5.5
pair_coeff 1 1 1 0.25 0 # A - A
pair_coeff 2 2 1 0.25 0 # U - U
pair_coeff 1 2 1 0.25 0 # A - U
pair_coeff 1 3 1 0.25 0 # A - M
pair_coeff 2 3 1 0.25 0 # U - M
pair_coeff 3 3 1 0.25 1 # M - M with 1/r^6 attraction (not Exponential)

# MD parameters
kspace_style ewald 1.0e-5
neighbor 2.0 bin
neigh_modify every 1 check yes
timestep 1.0 # 1 fs

