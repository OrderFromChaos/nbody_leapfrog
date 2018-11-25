"""
Convert to halos from different profiles
"""
import numpy as np
#from iccpy.cgs import msun, pc, G

c = 2.9979245800000e10 # speed of light, cm/s
pc = 3.0856775807000e18 # 1 parsec, cm
yr = 3.1556925200000e07 # year in s

m_p = 1.6726231100000e-24 # proton mass, g
G = 6.6725985000000e-08 # big G, cm^3 s^-2 g^-1
k_b = 1.3806581200000e-16 # boltzmann constant cm^2 s^-2 g K^-1 
 
sigma = 5.6705119000000e-5 # stefan boltzmann constant,  s^-3 g K^-4 
msun = 1.9889225000000e33 # solar mass, g  

h = 6.6260755400000E-27 # planck constant cm^2 s^-1 g 

kpc = 1e3 * pc
mpc = 1e3 * kpc

myr = 1e6 * yr

def solve_cc_from_rho0_rhoc(rho_ratio, cc_max=100.0):
    """
    rho_ratio - rho0/rho_crit where rho0 is the scale radius of the NFW density profile and
                rho_crit is the mean density of the universe (incl. dark energy)


    NFW density:
    rho(r) = rho0 / ( (r/rs) * (1+r/rs)^2)

    returns cc - concentration, i.e. r200/r_s such that the mean density is 200x the mean density

    """
    npts = 10000
    # r/r_s
    r_rs = (np.arange(npts)+1) / float(npts) * cc_max
    dens_ratio = 4*np.pi*rho_ratio* (np.log(1+r_rs) - r_rs/(1.0+r_rs)) / (4.0*np.pi * (r_rs**3)/3.0) 

    cc = np.interp(0.0, 200.0 - dens_ratio, r_rs) # search for 200*critical
    return cc

def halo_props(rho_s, r_s, H):
    """ solve for mass and concentration """
    # critical density
    rho_c = 3*(H*1e5/(1e6*pc))**2 / (8*np.pi *G) * ((1e3*pc)**3 / msun) # solar masses per kpc^3
    rho_ratio = rho_s/rho_c
    cc = solve_cc_from_rho0_rhoc(rho_ratio)

    m200 = 4*np.pi*rho_s*(r_s**3) * (np.log(1+cc) - cc/(1+cc)) 
    v200 = (m200 * msun * 10*G * H*1e5/(1e6*pc))**(1.0/3.0) / 1e5 # km/s M200 = v200^3  / (10 HG)  Creasey 2013
    r200 = ((m200*msun) * G ) / (v200*1e5)**2 / (1e3*pc) # r200 = G * m200 / v200^2 Mo Mao White

    return m200, v200, r200, cc

def rhos_rs_from_mass_conc_h(m200, cc, H):

    v200 = (m200 * msun * 10*G * H*1e5/(1e6*pc))**(1.0/3.0) / 1e5 # km/s M200 = v200^3  / (10 HG)  Creasey 2013
    r200 = ((m200*msun) * G ) / (v200*1e5)**2 / (1e3*pc) # r200 = G * m200 / v200^2 Mo Mao White
    
    r_s = r200/cc
    rho_c = 3*(H*1e5/(1e6*pc))**2 / (8*np.pi *G) * ((1e3*pc)**3 / msun) # solar masses per kpc^3
    
#    200.0 = 4*np.pi*(rho_s/rho_c)* (np.log(1+r_rs) - r_rs/(1.0+r_rs)) / (4.0*np.pi * (r_rs**3)/3.0) 
    rho_s = rho_c*200.0/(4*np.pi* (np.log(1+cc) - cc/(1.0+cc)) / (4.0*np.pi * (cc**3)/3.0))
    return rho_s, r_s
