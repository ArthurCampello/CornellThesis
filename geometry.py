# THIS FILE ENCODES GEOMETRIC PARAMETERS RELEVANT TO THE SINGLE-MODULE DETECTOR ASSEMBLY

###########
# imports #
###########

from inputs import * # import from inputs file
from math import ceil # import math ceil function
from math import floor # import math floor function

#######################################################################
# fixed geometry variables (hard-coded numbers are in units of 0.1mm) #
#######################################################################

scale= round(deltad/(0.1*10**(-3)),2) # calculates scaling factor from input deltad

mx = ceil(448/scale) # module x-dimension
my = ceil(178/scale) # module y-dimension
px = ceil(15/scale)  # excess cf-cf beyond module in x-direction
py = ceil(15/scale)  # excess cf-cf beyond module in y-direction
sz = ceil(6/scale)   # silicon z-dimension
gz = ceil(3/scale)   # glass z-dimension
cz = ceil(3/scale)   # carbon fiber z-dimension
fz = ceil(20/scale)  # carbon foam z-dimension
ts = ceil(240/scale) # spacing between tubes
tr = floor(10/scale) # radius of tubes
hsw = ceil(8/scale)  # module heat shunt width

sthicks = sum([ceil(i/scale) for i in thicks]) # calculates total thicknesses of input materials

###################################
# intermediate geometry variables #
###################################

# general dimensions
gridx = mx+2*px+2 # (module x-dimension) + 2*(x-padding) + 2*(1px shell)
gridy = my+2*py+2 # (module y-dimension) + 2*(y-padding) + 2*(1px shell)
gridz = sz+gz+cz+2*fz+sum(thicks)+2 #(carbin fiber z) + (glasss x) + (carbon foam z) + (sum of input layer thicknesses) + 2*(1px shell)

# carbon foam
foamx = [1,mx+2*px+1] # min and max x-values spanned by carbon foam
foamy = [1,my+2*py+1] # min and max y-values spanned by carbon foam
foamz = [1,2*fz+1] # min and max z-values spanned by carbon foam

# carbon fiber
fiberx = [1,mx+2*px+1] # min and max x-values spanned by carbon fiber
fibery = [1,my+2*py+1] # min and max y-values spanned by carbon fiber
fiberz = [foamz[1], foamz[1]+cz] # min and max z-values spanned by carbon fiber

# input materials
inx = [1+px,1+px+mx] # min and max x-values spanned by input layers
iny = [1+py,1+py+my] # min and max y-values spanned by input layers
inz = [fiberz[1],fiberz[1]+sthicks] # min z-value spanned by input layrs

# glass
glassx = [1+px,1+px+mx] # min and max x-values spanned by glass
glassy = [1+py,1+py+my] # min and max y-values spanned by glass
glassz = [inz[1], inz[1] + gz] # min and max z-values spanned by glass

# silicon
six = [1+px,1+px+mx] # min and max x-values spanned by silicon
siy = [1+py,1+py+my] # min and max y-values spanned by silicon
siz = [glassz[1], glassz[1] + sz] # min and max z-values spanned by silicon

# tubing parameters
ctb1x = 1+px+int(mx/2)-int(ts/2) # x-position of center of tube 1
ctb2x = 1+px+int(mx/2)+int(ts/2) # x-position of center of tube 2
ctbz = 1+fz # z-position of center of both tubes

# module heat shunt parameters
vmod = (six[1]-six[0])*(siy[1]-siy[0])*(siz[1]-siz[0]) # volume of the module (in pixels)
vshunt = (six[1]-six[0])*(hsw)*(siz[1]-siz[0]) # volume of the shunts (in pixels)

# variables for Numba compatibility (these variables are needed for Numba not to throw an error)
foamz0 = foamz[0]
foamz1 = foamz[1]
fiberz0 = fiberz[0]
fiberz1 = fiberz[1]
inx0 = inx[0]
inx1 = inx[1]
iny0 = iny[0]
iny1 = iny[1]
inz0 = inz[0]
inz1 = inz[1]
glassz0 = glassz[0]
glassz1 = glassz[1]
siz0 = siz[0]
siz1 = siz[1]
