#!/usr/bin/env python
import numpy as np

SB = 5.6696E-8 # Stefan-Boltzmann [W/m^2/K^4]
S_0 = 1368 # Solar constant [W/m^2]
R_E = 6371E3 # Radius of Earth [m]
A_E = 4*np.pi*R_E**2 # Total surface area of Earth [m^2]

emissivity = 1 # Total emissivity (play around with this?)
transmissivity = 0.63 # Atmospheric transmissivity
alpha_sky = 0.2 # Atmospheric albedo

alpha_land = 0.4 # Albedo of land surface
alpha_ocean = 0.1 # Albedo of ocean surface
alpha_ice = 0.6 # Albedo of ice

rho_land = 2500 # Density of land surface [kg/m^3]
rho_ocean = 1028 # Density of ocean surface [kg/m^3]
rho_ice = 900 # Density of ice [kg/m^3]

Z_land = 1.0 # Thermal scale depth for land [m]
Z_ocean = 70.0 # Thermal scale depth for ocean [m]
Z_ice = 1.0 # Thermal scale depth for ice [m]

c_land = 790 # Specific heat capacity for land [J/kg/K]
c_ocean = 4187 # Specific heat capacity for water [J/kg/K]
c_ice = 2060 # Specific heat capacity for ice [J/kg/K]

# Vector of products of density, thermal length scale, and specific heat capacity for land, ocean, ice
LOI_rho_Z_c = np.array([rho_land*Z_land*c_land, rho_ocean*Z_ocean*c_ocean, rho_ice*Z_ice*c_ice])
# Vector of albedos for land, ocean, ice
LOI_alpha = np.array([alpha_land, alpha_ocean, alpha_ice])

# Land, ocean, ice fractions for each zone as a vector. Columns are for [land, ocean, ice], rows correspond to each zone
fractions = np.genfromtxt('data/model-zonal-fractions.csv', delimiter=',')

# Vector constants dependent on zone
k = np.genfromtxt('data/thermalexchangecoeffs.csv', delimiter=',') # Thermal exchange coefficient [W/m/K]
a = np.zeros(8)
a[1:7] = 0.5*np.array(list(map(lambda i: np.sin(np.pi*(-0.5+(i+1)/6)) - np.sin(np.pi*(-0.5+i/6)), range(6)))) # Zonal area fractions vector
gamma = np.array([0, 0.1076, 0.2277, 0.3045, 0.3045, 0.2277, 0.1076, 0]) # Zonal geometric fractions vector
L = np.zeros(7)
L[1:6] = 2*np.pi*R_E*np.array(list(map(lambda i: np.cos(np.pi*(-0.5+i/6)), range(1,6)))) # Boundary lengths vector

# Area averaged properties for each zone as vectors
rho_c_Z = np.pad(np.dot(fractions, LOI_rho_Z_c), (1, 1), constant_values=(1, 1)) # Area averaged product of density, specific heat capacity, thermal scale depth 
alpha = np.pad(np.dot(fractions, LOI_alpha), (1, 1)) # Area averaged albedo