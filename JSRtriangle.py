"""
Created on Fri Aug 17 2018
@author: RodgerCornell
"""
import thermo as therm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# JSR OPERATIONAL SPACE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#------------------------------------------------------------------------------
# Dimensions of the Reactor
#------------------------------------------------------------------------------
#Reactor = 0.029   # Burke Lab
#Nozzle = 300      # Burke Lab 
Reactor = 0.070    # Test (Matras R1)
Nozzle = 1000      # Test (Matras R1)

# Reactor Temperature [K]
Temp =  20 + 273.15

# Residence Time [s]
tau = 5.0

# Adimensional Constant
# *(we need to find more data on this constant. thers is supposedly a paper
# that provides A as a function of temperature) 
A = 0.3 #np.pi/4

#------------------------------------------------------------------------------
# Species Information/Thermo Data
#------------------------------------------------------------------------------
# Stoichiometric Air/Methane (example)
#  (1)CH4 +  (2)(02 + 3.76N2) = (2)H2O + (1)CO2 + (7.52)N2
#------------------------------------------------------------------------------
N2 = therm.Chemical('n2',T=Temp)            # Test (Matras R1)
O2 = therm.Chemical('o2',T=Temp)            # Test (Matras R1)
Ar = therm.Chemical('argon',T=Temp)         # Test (Matras R1)

# (N2, O2, Ar)
# Mole Fractions []
n = np.array([0.7809,0.2095,0.0096])

# Molecular Mass [g/mol]
M = np.array([N2.MW,O2.MW,Ar.MW])

# Mass Fractions []
# *(mole fraction below if given mass fraction) 
y = (n*M)/(sum(n*M))
#n = (y/M)/sum(y/M)

# Universal Gas Constant [J/(K*mol)]
R_uni = 8.3144598 

# Cp [J/(K*g)]
Cp = np.array([N2.Cpg/1000,O2.Cpg/1000,Ar.Cpg/1000])

# Mass Density [g/m^3]
rho = np.array([N2.rhog*1000,O2.rhog*1000,Ar.rhog*1000]) 

# Dynamic Viscosity [(g/(m*s^2))*s]
nu = np.array([N2.mug*1000,O2.mug*1000,Ar.mug*1000])

#------------------------------------------------------------------------------
# Mixture Information/Thermo Data
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

#------------------------------------------------------------------------------
# The Three Limits
#------------------------------------------------------------------------------
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

# RECYCLING RATE
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
    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOTTING THE RESULTS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#------------------------------------------------------------------------------
# Finding the Corners of the Triangle to Bound an Operational Space
#------------------------------------------------------------------------------
point1d = (920*nu_mix)/(3*rho_mix*A*c_mix)
point1R = ((3*tau*(point1d**2)*c_mix)/4)**(1/3)

point2d = (3*tau*c_mix*((np.pi)**3)*(A**3))/864000
point2R = (60*point2d)/(np.pi*A)

point3d = ((230*tau*nu_mix*(np.pi**3)*(A**2))/(rho_mix*(60**3)))**(1/2)
point3R = ((230*tau*nu_mix*point3d)/(rho_mix*A))**(1/3)

#------------------------------------------------------------------------------
# Setting Plot Boundaries 
#------------------------------------------------------------------------------
plt.xlim([point1R - 0.3*point1R,point2R + 0.3*point2R])                                    
plt.ylim([point1d - 0.3*point1d,point2d + 0.3*point2d])                                

#------------------------------------------------------------------------------
# Setting the aspect ratio of the plot
#------------------------------------------------------------------------------
plt.gca().set_aspect(aspect=0.5,adjustable='box')

#------------------------------------------------------------------------------
# Combining Turbulence and Recycling lines to make shading possible 
#------------------------------------------------------------------------------
line4 = np.maximum(Rt,Rr)
# Bounding this new line using the triangle corners
line4 = line4[np.logical_and(line4>point1R,line4<point2R)]

#------------------------------------------------------------------------------
# Shading the triangle and locating the plot text
#------------------------------------------------------------------------------
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
plt.fill_betweenx(step,line4,Rs_cut,color='grey',alpha='0.2')     

# Plotting the data on a log-log scale     
plt.loglog(Reactor,Nozzle*10**(-6),"o",color='k')
plt.loglog(Rt,dt,color='b', linestyle='--')
plt.loglog(Rr,dr,color='r', linestyle='-.')
plt.loglog(Rs,ds,color='g', linestyle='-')
    
#------------------------------------------------------------------------------
plt.title('JSR Operational Limits', color='k')
plt.ylabel('Nozzle Diameter [m]', color='k')
plt.xlabel('Reactor Radius [m]', color='k')
plt.legend(['Reactor','Turbulence Limit','Recycling Limit','Sonic Limit'],
           fancybox=True,loc='lower right')
# Adding tick marks to all sides of the plot 
plt.tick_params(axis='both',which='both',direction='in',width=1.25,top=True,
                right=True)
plt.savefig('JSRplot.jpg',dpi=1200,bbox_inches='tight')

#------------------------------------------------------------------------------


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MATRAS VALIDATION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#------------------------------------------------------------------------------
# Pulling in Digitized Matras Data for Validation
#------------------------------------------------------------------------------
#df1 = pd.read_excel("Matras_1s.xlsx")
#plt.loglog(df1['Rt [m]'],df1['dt [m]'],'x',color='b')
#plt.loglog(df1['Rr [m]'],df1['dr [m]'],'x',color='r')
#plt.loglog(df1['Rs [m]'],df1['ds [m]'],'x',color='g')
#plt.legend(['Turbulence','Recycling','Sonic','Matras_T','Matras_R','Matras_S']
#           ,fancybox=True,loc='lower right')
#plt.savefig('Matras_1s_validation.jpg',dpi=1200,bbox_inches='tight')

#df5 = pd.read_excel("Matras_5s.xlsx")
#plt.loglog(df5['Rt [m]'],df5['dt [m]'],'x',color='b')
#plt.loglog(df5['Rr [m]'],df5['dr [m]'],'x',color='r')
#plt.loglog(df5['Rs [m]'],df5['ds [m]'],'x',color='g')
#plt.legend(['Turbulence','Recycling','Sonic','Matras_T','Matras_R','Matras_S']
#           ,fancybox=True,loc='lower right')
#plt.savefig('Matras_5s_validation.jpg',dpi=1200,bbox_inches='tight')





 