import numpy as np

# temperature and heat
at = -10  # initial ambient temperature in degrees C
tt = -30 # tube temperature in degrees C

deltad = 0.3*10**(-3)
foamtype = 2 # Can be either 1, 2, or 3 with 1 corresponding to 0.09g/cc, 2 to 0.19g/cc, 3 to 0.35g/cc

mpsr = 0.6 # module power shunt ratio 

#thicks = [1,1,1] # thicknesses of intermediate materials in units of 0.1mm
#alphas = [np.array([0.00004]),np.array([.000072507]),np.array([0.00004])] # alphas of intermediate materials in units of 0.1mm

thicks = [1] # thicknesses of intermediate materials in units of 0.1mm
alphas = [np.array([0.00008])] # alphas of intermediate materials in units of 0.1mm

#thicks = [] # thicknesses of intermediate materials in units of 0.1mm
#alphas = [np.array([0])] # alphas of intermediate materials in units of 0.1mm

# this appears as a list of lists (e.g. [np.array([1]),np.array([6]),np.array([2,5])]) where each 
# sublist is of size one for materials of isotropic thermal conductivities and
# size 3 for materials with anisotropic thermal conductivities
# [[0 by default]]