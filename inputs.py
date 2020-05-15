###########
# imports #
###########

import numpy as np # import numpy

# temperature and heat
at = -10 # initial ambient temperature in degrees C (default = -10)
tt = -30 # tube temperature in degrees C (default = -30)

deltad = 0.3*10**(-3) # the pixel width (default = 0.3*10**(-3))
foamtype = 2 # Can be either 1, 2, or 3 with 1 corresponding to 0.09g/cc, 2 to 0.19g/cc, 3 to 0.35g/cc

mpsr = 0.6 # module power shunt ratio (default = 0.6, shown in Figure 4.1 a) in thesis)

alpha_contact = 0.00002 # thermal diccusivity of one-pixel thick 'contact' layer representing adhesive or orginary contact (default =0.00002)

# inputs related to materials between carbon fiber and module
thicks = [] # thicknesses of intermediate materials in units of 0.1mm
alphas = [np.array([0])] # alphas of intermediate materials in units of 0.1mm (placed in nested arrays for thermal anisotropy possibility)

# inputs corresponding to one contact layer between carbon fiber and module (Assembly 1 in Figure 5.1 in thesis)
#thicks = [1] # thicknesses of intermediate materials in units of 0.1mm
#alphas = [np.array([alpha_contact])] # alphas of intermediate materials in units of 0.1mm

# inputs corresponding to two contact layers between carbon fiber and module with AlN layer in between (Assembly 3 in Figure 5.1 in thesis)
#thicks = [1,1,1] # thicknesses of intermediate materials in units of 0.1mm
#alphas = [np.array([alpha_contact]),np.array([.000072507]),np.array([alpha_contact])] # alphas of intermediate materials in units of 0.1mm