# units μm, hr, pg 

# gravity constant 
g_SI = 9.81 # m/s^2
g = g_SI * 1e6 * 60^4 # μm/hr^2
g_per_second = g_SI * 1e6 # μm/s^2

# densities 
c_water = 1 # pg/μm^3
c_cell = 1.05 # pg/μm^3

# mass of a cell 
r = 5 # μm
m = 4/3 * π * r^3 * (c_cell-c_water) # pg

# stokes law, with dynamic viscosity of water * 100 times 
# see, The Physics of Living Systems, Page 239:
μ_SI = 1.0 * 100 # mP s = kg / (m s²)    
μ = μ_SI * 1e12 * 1e-6 * 60^4 # pg/(μm hr²)
η = 6 * π * μ * r # μm^2/hr


# speed 
v = g * m / η # μm/hr

v_per_minute = v / 60 # μm/min

