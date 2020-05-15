# THIS FILE ENCODES MATERIAL PARAMETERS RELEVANT TO THE SINGLE-MODULE DETECTOR ASSEMBLY

###########
# imports #
###########

from inputs import * # import from inputs file
import numpy as np # import numpy

###########################################
# fixed material variables (with sources) #
###########################################

# carbon foam
if foamtype == 1:
    dens_cfo = 90 #(kg/m^3) # according to internal MaterialData.xml
    k_cfo = 10 #(W/(m*K)) # according to internal MaterialData.xml
elif foamtype == 2:
    dens_cfo = 190 #(kg/m^3) # according to internal MaterialData.xml
    k_cfo = 29 #(W/(m*K)) # according to internal MaterialData.xml
elif foamtype == 3:
    dens_cfo = 350 #(kg/m^3) # according to internal MaterialData.xml
    k_cfo = 75 #(W/(m*K)) # according to internal MaterialData.xml
else:
    print('Error: Invalid Foam Type')
cp_cfo = 180 #(J/(K*kg))
alpha_cfo = (k_cfo)/(cp_cfo*dens_cfo) # calculate foam thermal diffusivity

# carbon fiber and epoxy (AS4/3501-6 Composite) 
cp_cfi = 887 #(J/(K*kg)) 
dens_cfi = 1578 #(kg/m^3) # according to internal MaterialData.xml
k_cfi = np.array([15.7, 15.7, 0.687]) #(W/(m*K))
alpha_cfi = k_cfi/(cp_cfi*dens_cfi) # calculate fiber thermal diffusivity

# glass
cp_glass = 840 #(J/(K*kg)) https://www.engineeringtoolbox.com/specific-heat-solids-d_154.html
dens_glass = 2600 #(kg.m^3) https://www.engineeringtoolbox.com/density-solids-d_1265.html
k_glass = 1.01 #(W/(m*K)) # https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
def k_glass_T(T):
    return 1.0143+T*3.57*10**(-4)
alpha_glass = (k_glass)/(cp_glass*dens_glass) # calculate glass thermal diffusivity

# silicon
cp_si = 705 #(J/(K*kg)) https://www.engineeringtoolbox.com/specific-heat-solids-d_154.html
dens_si = 2336 #(kg/m^3) # according to internal MaterialData.xml
k_si = 148 #(W/(m*K)) # according to internal MaterialData.xml
def k_si_T(T):
    return 165.406-T*0.625
alpha_si = (k_si)/(cp_si*dens_si) # calculate silicon thermal diffusivity

# aluminum nitride
cp_aln = 750 #(J/(K*kg)) https://precision-ceramics.com/materials/aluminum-nitride/
dens_aln = 3310 #(kg/m^3) https://precision-ceramics.com/materials/aluminum-nitride/
k_aln = 180 #(W/(m*K)) https://precision-ceramics.com/materials/aluminum-nitride/
alpha_aln = k_aln/(cp_aln*dens_aln) # calculate aluminum nitride thermal diffusivity

##################################################
# intermediate material and simulation variables #
##################################################

tmaxr = 1/6 # target maximum r (can be up to 1/6)
maxalphasins = max([max(i) for i in alphas]) # find the maximum thermal diffusivity of materials in input layers
dt = tmaxr/(max(alpha_cfo, max(alpha_cfi), alpha_glass, alpha_si, maxalphasins, alpha_aln)/deltad**2) # calculate time scale

r_cfo = dt*alpha_cfo/deltad**2 # calculate r-parameter of carbon foam

r_cfi = dt*alpha_cfi/deltad**2 # calculate r-parameter of carbon fiber

def r_glass_T(T): # calculate r-parameter of glass as a function of temperature
    return dt*((k_glass_T(T))/(cp_glass*dens_glass))/deltad**2

def r_si_T(T): # calculate r-parameter of silicon as a function of temperature
    return dt*((k_si_T(T))/(cp_si*dens_si))/deltad**2

r_aln = dt*alpha_aln/deltad**2 # calculate r-parameter of aluminum nitride

# calculates all r-parameters of input materials
r_ins = []
for i in alphas:
    r_ins.append(dt*i/deltad**2)
