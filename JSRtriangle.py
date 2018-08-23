# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 08:55:49 2018

@author: RodgerCornell
"""
import thermo as therm
import matplotlib.pyplot as plt
import numpy as np
#------------------------------------------------------------------------------
# JSR OPERATIONAL SPACE
#------------------------------------------------------------------------------
Reactor = 0.020 # Test (Van Geem)    0.029 #[m]
Nozzle = 900    # Test (Van Geem)    300    #[um]

# Reactor Temperature [K]
Temp =  900 # Test (Van Geem)      800
# Residence Time [s]
tau = 0.2
# Adimensional Constant
A = np.pi/4 #0.3

#------------------------------------------------------------------------------
# SPECIES INFO
#------------------------------------------------------------------------------
# Stoichiometric Air/Methane
#  (1)CH4 +  (2)(02 + 3.76N2) = (2)H2O + (1)CO2 + (7.52)N2
#------------------------------------------------------------------------------
C2H6 = therm.Chemical('ethane',T=Temp)                      # Test (Van Geem)
H2O = therm.Chemical('water',T=Temp)                        # Test (Van Geem)

y = ([0.769,0.231])                                         # Test (Van Geem)

CH4 = therm.Chemical('methane',T=Temp)
O2 = therm.Chemical('oxygen',T=Temp)
N2 = therm.Chemical('nitrogen',T=Temp)

# (CH4, O2, N2)
## Mole Fractions []
#n = np.array([0.095057,0.190114,0.714829])

# Molecular Mass [g/mol]
M = np.array([C2H6.MW,H2O.MW])                              # Test (Van Geem)
#M = np.array([CH4.MW,O2.MW,N2.MW])

## Mass Fractions []
#y = (n*M)/(sum(n*M))
n = (y/M)/sum(y/M)                                          # Test (Van Geem)

# Universal Gas Constant [J/(K*mol)]
R_uni = 8.3144598 

# Cp [J/(K*g)]
Cp = np.array([C2H6.Cpg/1000,H2O.Cpg/1000])                 # Test (Van Geem)
#Cp = np.array([CH4.Cpg/1000,O2.Cpg/1000,N2.Cpg/1000])

# Mass Density [g/m^3]
rho = np.array([C2H6.rhog*1000,H2O.rhog*1000])              # Test (Van Geem) 
#rho = np.array([CH4.rhog*1000,O2.rhog*1000,N2.rhog*1000]) 

# Dynamic Viscosity [(g/(m*s^2))*s]
nu = np.array([C2H6.mug*1000,H2O.mug*1000])                 # Test (Van Geem)
#nu = np.array([CH4.mug*1000,O2.mug*1000,N2.mug*1000])

#------------------------------------------------------------------------------
# MIXTURE INFO
#------------------------------------------------------------------------------
# Molecular Mass [g/mol]
M_mix = sum(n*M)
# Gas Constant [J/(K*g)]
R_mix = (R_uni)/M_mix
# Mass Density [g/m^3]
rho_mix = sum(n*rho)
# Dynamic Viscosity [g/(m*s)]
nu_mix = sum(n*nu)
# Speed of Sound (T,P) [m/s]
k_mix = (sum(y*Cp))/((sum(y*Cp)-R_mix))
c_mix = np.sqrt(k_mix*R_mix*Temp*1000)

# TURBULENCE
# The jets flowing from the nozzles must be turbulent. This imposes an upper
# limit of the residence time. Putting this in terms of Reactor Radius (R)
# and Nozzle Diameter (d) and making it an equality...
dt=[]
Rt=[]
for x in range(1,Nozzle*1000,1):
    x = x*10**(-6)
    dt.append(x)
    Rt.append(((230*tau*nu_mix*x)/(rho_mix*A))**(1/3))

# RECYCLING Rate
# This is the ratio of residence time to the time it takes to carry all gas in 
# the reactor from one nozzle to the adjacent one. Experimentally observed that
# recycling rate must be >30 for good gas-phase mixing. This creates a lower
# limit on the ratio of Reactor Radius to Nozzle Diameter. Writing in terms of 
# Reactor Radius (R) and Nozzle Diameter (d) and making it an equality...   
dr=[]
Rr=[]
for x in range(1,Nozzle*1000,1):
    x = x*10**(-6)
    dr.append(x)
    Rr.append((60*x)/(np.pi*A))
    
# SONIC LIMIT
# Velocity of the gas exiting the nozzle must be less than the speed of sound
# at the T and P of the reaction. This will create a lower limit on the 
# residence time. Putting this in terms of Reactor Radius (R) and Nozzle
ds=[]
Rs=[]
for x in range(1,Nozzle*1000,1):
    x = x*10**(-6)
    ds.append(x)    
    Rs.append(((3*tau*(x**2)*c_mix)/4)**(1/3))
    
#------------------------------------------------------------------------------
# PLOTTING RESULTS
#------------------------------------------------------------------------------
# Finding the ends of the triangle to bound the shaded region
point1d = (920*nu_mix)/(3*rho_mix*A*c_mix)
point1R = ((3*tau*(point1d**2)*c_mix)/4)**(1/3)

point2d = (3*tau*c_mix*((np.pi)**3)*(A**3))/864000
point2R = (60*point2d)/(np.pi*A)

point3d = ((230*tau*nu_mix*(np.pi**3)*(A**2))/(rho_mix*(60**3)))**(1/2)
point3R = ((230*tau*nu_mix*point3d)/(rho_mix*A))**(1/3)

#plt.loglog(point1R,point1d,"o",color='y')
#plt.loglog(point2R,point2d,"o",color='y')
#plt.loglog(point3R,point3d,"o",color='y')

#------------------------------------------------------------------------------
plt.xlim([point1R - 0.2*point1R,point2R + 0.2*point2R])                                    
plt.ylim([point1d - 0.2*point1d,point2d + 0.2*point2d])                                

#------------------------------------------------------------------------------
# Combining two lines into one to make shading possible with 'fill_between' 
line4 = np.maximum(Rt,Rr)
# Cutting down the line
line4 = line4[np.logical_and(line4>point1R,line4<point2R)]
#line4 = line4[np.logical_and(x>point1R,x<point2R)]
#------------------------------------------------------------------------------
# Using the triangle ends to limit the shading and locate the plot text
step=[]
ds_cut=[]
Rs_cut=[]
dr_cut=[]
Rr_cut=[]
dt_cut=[]
Rt_cut=[]
d=0
for x in range(1,Nozzle*1000,1):
    y = x
    x = x*10**(-6)
    if point1d < x < point2d:
        step.append(x)
        ds_cut.append(ds[y])
        Rs_cut.append(Rs[y])
    if point3d < x < point2d:
        dr_cut.append(dr[y])
        Rr_cut.append(Rr[y])
    if point1d < x < point3d:
        dt_cut.append(dt[y])
        Rt_cut.append(Rt[y])

# Finding logarithmic means to plot labels
Mean13x = (point3R-point1R)/(np.log(point3R)-np.log(point1R))
Mean13y = (point3d-point1d)/(np.log(point3d)-np.log(point1d))  
Mean32x = (point3R-point2R)/(np.log2(point3R)-np.log2(point2R))
Mean32y = (point3d-point2d)/(np.log2(point3d)-np.log2(point2d))  
Mean12x = (point2R-point1R)/(np.log(point2R)-np.log(point1R))
Mean12y = (point2d-point1d)/(np.log2(point2d)-np.log2(point1d)) 
     
plt.loglog(Reactor,Nozzle*10**(-6),"o",color='k')
plt.loglog(Rt,dt,color='c', linestyle='--')
plt.loglog(Rr,dr,color='r', linestyle='-.')
plt.loglog(Rs,ds,color='g', linestyle='-')
plt.legend(['Reactor', 'Turbulence','Recycling','Sonic'],fancybox=True,loc='lower right')
        
plt.fill_betweenx(step,line4,Rs_cut,color='grey',alpha='0.2')      
#------------------------------------------------------------------------------
plt.title('JSR Operational Limits', color='k')
plt.ylabel('Nozzle Diameter (m)', color='k')
plt.xlabel('Reactor Radius (m)', color='k')

plt.tick_params(axis='both',which='both',direction='in',width=1.25,top=True,right=True)

#plt.text(Mean12x,Mean12y,'Sonic',color='g')
#plt.text(Mean13x,Mean13y,'Turbulence',color='b')
#plt.text(Mean32x,Mean32y,'Recycling',color='r')

plt.savefig('whatever.jpg',dpi=1200,bbox_inches='tight')
#------------------------------------------------------------------------------ 