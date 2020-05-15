###########
# imports #
###########

from inputs import *
import numpy as np

###########################################
# fixed material variables (with sources) #
###########################################

# carbon foam (K9 CVD Foam 0.19g/cc)
if foamtype == 1:
    dens_cfo = 90 #(kg/m^3) # according to MaterialData.xml
    k_cfo = 10 #(W/(m*K)) # according to MaterialData.xml
elif foamtype == 2:
    dens_cfo = 190 #(kg/m^3) # according to MaterialData.xml
    k_cfo = 29 #(W/(m*K)) # according to MaterialData.xml
elif foamtype == 3:
    dens_cfo = 350 #(kg/m^3) # according to MaterialData.xml
    k_cfo = 75 #(W/(m*K)) # according to MaterialData.xml
else:
    print('Error: Invalid Foam Type')
cp_cfo = 180 #(J/(K*kg)) # use file:///C:/Users/campe/Downloads/dayton1271365522.pdf as source
alpha_cfo = (k_cfo)/(cp_cfo*dens_cfo)

# carbon fiber and epoxy (AS4/3501-6 Composite) file:///C:/Users/campe/Downloads/ohiou1228761926%20(1).pdf 
cp_cfi = 887 #(J/(K*kg)) 
dens_cfi = 1578 #(kg/m^3) # according to MaterialData.xml
k_cfi = np.array([15.7, 15.7, 0.687])
alpha_cfi = k_cfi/(cp_cfi*dens_cfi)

# glass
cp_glass = 840 #(J/(K*kg)) https://www.engineeringtoolbox.com/specific-heat-solids-d_154.html
dens_glass = 2600 #(kg.m^3) https://www.engineeringtoolbox.com/density-solids-d_1265.html
k_glass = 1.01 #(W/(m*K)) # https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
def k_glass_T(T):
    return 1.0143+T*3.57*10**(-4)
alpha_glass = (k_glass)/(cp_glass*dens_glass) # should be about 8*10^(-5) according to http://www.ioffe.ru/SVA/NSM/Semicond/Si/thermal.html

# silicon
cp_si = 705 #(J/(K*kg)) https://www.engineeringtoolbox.com/specific-heat-solids-d_154.html
dens_si = 2336 #(kg/m^3) # according to MaterialData.xml
k_si = 148 #(W/(m*K)) # according to MaterialData.xml
def k_si_T(T):
    return 165.406-T*0.625
alpha_si = (k_si)/(cp_si*dens_si) # should be about 8*10^(-5) according to http://www.ioffe.ru/SVA/NSM/Semicond/Si/thermal.html

# aluminum nitride
cp_aln = 750 #(J/(K*kg)) 
dens_aln = 3310 #(kg/m^3) # according to MaterialData.xml
k_aln = 180
alpha_aln = k_aln/(cp_aln*dens_aln)

##################################################
# intermediate material and simulation variables #
##################################################

tmaxr = 1/6 # target maximum r (can be up to 1/6)
maxalphasins = max([max(i) for i in alphas])
dt = tmaxr/(max(alpha_cfo, max(alpha_cfi), alpha_glass, alpha_si, maxalphasins, alpha_aln)/deltad**2)

r_cfo = dt*alpha_cfo/deltad**2
r_cfi = dt*alpha_cfi/deltad**2

def r_glass_T(T):
    return dt*((k_glass_T(T))/(cp_glass*dens_glass))/deltad**2

def r_si_T(T):
    return dt*((k_si_T(T))/(cp_si*dens_si))/deltad**2

r_aln = dt*alpha_aln/deltad**2

r_ins = []
for i in alphas:
    r_ins.append(dt*i/deltad**2)