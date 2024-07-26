#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for quick validation of results to test e.g. single precision or optimizations. 
Compare the fluxes and heating rates between two ecRad output files e.g.:
./compare_flux_hr.py ecrad_meridian_ecckd_spartacus_out.nc ecrad_meridian_ecckd_spartacus_out_REFERENCE.nc 

@author: Peter Ukkonen
"""


import os
import sys
import numpy as np
from netCDF4 import Dataset


def mae(predictions,targets):
    diff = predictions - targets
    return np.mean(np.abs(diff))

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def bias(predictions, targets):
    return np.mean(predictions-targets)

def calc_heatingrate(fluxup,fluxdn,p):
    F = fluxdn - fluxup
    dF = F[:,1:] - F[:,0:-1] 
    dp = p[:,1:] - p[:,0:-1] 
    dFdp = dF/dp
    g = 9.81 # m s-2
    cp = 1004 # J K-1  kg-1
    dTdt = -(g/cp)*(dFdp) # K / s
    dTdt_day = (24*3600)*dTdt
    return dTdt_day


if len(sys.argv) > 3:
    print('You have specified too many arguments')
    sys.exit()

if len(sys.argv) < 3:
    print('You need two files to compare')
    sys.exit()

fpath           = sys.argv[1]
fpath_ref       = sys.argv[2]

# fpath_ref   = os.path.splitext(fpath)[0] + '_REFERENCE.nc'

if not os.path.isfile(fpath):
    print('The file specified does not exist ({})'.format(fpath))
    sys.exit()

if not os.path.isfile(fpath_ref):
    print('Corresponding reference file ({}) not found'.format(fpath_ref))
    sys.exit()
    

dat             =  Dataset(fpath)
dat_ref         =  Dataset(fpath_ref)

pres   = dat_ref.variables['pressure_hl'][:,:].data
pres_lay = 0.5 * (pres[:,1:] + pres[:,0:-1])

ncol = pres.shape[0]
nlev = pres.shape[1]

rlu_ref = dat_ref.variables['flux_up_lw'][:,:].data
rld_ref = dat_ref.variables['flux_dn_lw'][:,:].data
rsu_ref = dat_ref.variables['flux_up_sw'][:,:].data
rsd_ref = dat_ref.variables['flux_dn_sw'][:,:].data


if 'flux_up_lw' in dat.variables:
    rlu     = dat.variables['flux_up_lw'][:,:].data
    rld     = dat.variables['flux_dn_lw'][:,:].data
    do_lw = True
else:
    do_lw = False
if 'flux_up_sw' in dat.variables:
    rsu     = dat.variables['flux_up_sw'][:,:].data
    rsd     = dat.variables['flux_dn_sw'][:,:].data
    do_sw = True
else:
    do_sw = False
    
if 'flux_dn_sw_clear' in dat.variables:
    rsu_clear     = dat.variables['flux_up_sw_clear'][:,:].data
    rsd_clear     = dat.variables['flux_dn_sw_clear'][:,:].data
    rsu_clear_ref  = dat_ref.variables['flux_up_sw_clear'][:,:].data
    rsd_clear_ref  = dat_ref.variables['flux_dn_sw_clear'][:,:].data
    do_sw_clear = True
    
else:
    do_sw_clear = False
    
if 'flux_dn_lw_clear' in dat.variables:
    rlu_clear     = dat.variables['flux_up_lw_clear'][:,:].data
    rld_clear     = dat.variables['flux_dn_lw_clear'][:,:].data
    rlu_clear_ref     = dat_ref.variables['flux_up_lw_clear'][:,:].data
    rld_clear_ref     = dat_ref.variables['flux_dn_lw_clear'][:,:].data
    do_lw_clear = True
else:
    do_lw_clear = False
    
dTdt_lw_ref     = calc_heatingrate(rlu_ref,rld_ref,pres)
dTdt_sw_ref     = calc_heatingrate(rsu_ref,rsd_ref,pres)
dTdt_lw         = calc_heatingrate(rlu,rld,pres)
dTdt_sw         = calc_heatingrate(rsu,rsd,pres)

dat.close()
dat_ref.close()

# print("shape: {}".format(rld.shape))
print ("-- Comparing the fluxes and heating rates in the two files")
print ("-- one of the files may be a reference computation (calling differences 'error')")

if do_lw:
    
    print("--------- COMPARING LONGWAVE (LW) OUTPUTS -----------------")
    
    # Net flux = downwelling - upwelling
    flux_lw_ref = rld_ref - rlu_ref
    flux_lw     = rld - rlu
    
    diff = np.max(np.abs(flux_lw - flux_lw_ref))
    print("Max. absolute error in LW net fluxes:")
    
    if (diff < 0.3):
          print("{:0.3f} = less than 0.3 W/m2 -------> Check PASSED!".format(diff))
    else:
          print("{:0.3f} = more than 0.3 W/m2 -------> Check NOT PASSED".format(diff))
    
    
    diff = mae(flux_lw, flux_lw_ref)
    print("Mean absolute error in LW net fluxes:")
    
    if (diff < 0.05):
          print("{:0.3f} = less than 0.05 W/m2 ------> Check PASSED!".format(diff))
    else:
          print("{:0.3f} = more than 0.05 W/m2 ------> Check NOT PASSED".format(diff))
    
    
    diff =  mae(rlu, rlu_ref)
    print("Mean absolute error in upwelling LW flux at top-of-atmosphere:")
    if (diff < 0.05):
         print("{:0.3f} = less than 0.05 W/m2 ------> Check PASSED!".format(diff))
    else:
         print("{:0.3f} = more than 0.05 W/m2 ------> Check NOT PASSED".format(diff))
    
    if do_lw_clear:
        flux_lw_clear_ref = rld_clear_ref - rlu_clear_ref
        flux_lw_clear     = rld_clear - rlu_clear
    
        diff = np.max(np.abs(flux_lw_clear-flux_lw_clear_ref))
    
        print("Max error in clear-sky longwave net fluxes:")
        print("{:0.3f}".format(diff))
        
    print("----- HEATING RATES ------- ")
    
    diff = bias(dTdt_lw, dTdt_lw_ref)
    print("Bias in LW heating rate:")
    print("{:0.4f}".format(diff))

    inds_trop = pres_lay > 10000.0
    inds_strat = pres_lay < 10000.0

    diff = bias(dTdt_lw[inds_trop], dTdt_lw_ref[inds_trop])
    print("Bias in LW heating rate (troposphere):")
    print("{:0.4f}".format(diff))
	    
    diff = mae(dTdt_lw, dTdt_lw_ref)
    print("Mean absolute error in LW heating rate:")
    
    if (diff < 0.02):
          print("{:0.3f} = less than 0.02 K/day ------> Check PASSED!".format(diff))
    else:
          print("{:0.3f} = more than 0.02 K/day ------> Check NOT PASSED".format(diff))
    

    diff = np.max(np.abs(dTdt_lw-dTdt_lw_ref))
    print("Max. absolute error in LW heating rates:")
    if (diff < 0.5):
         print("{:0.3f} = less than 0.5 K/day ------> Check PASSED!".format(diff))
    else:
         print("{:0.3f} = more than 0.5 K/day ------> Check NOT PASSED".format(diff))

    inds_trop = pres_lay > 10000.0
    inds_strat = pres_lay < 10000.0
    diff = np.max(np.abs(dTdt_lw[inds_trop]-dTdt_lw_ref[inds_trop]))
    print("Max. absolute error in LW heating rates (troposphere):")
    if (diff < 0.2):
         print("{:0.3f} = less than 0.2 K/day ------> Check PASSED!".format(diff))
    else:
         print("{:0.3f} = more than 0.2 K/day ------> Check NOT PASSED".format(diff))



if do_sw:
    print("--------- COMPARING SHORTWAVE (SW) OUTPUTS -----------------")
    
    
    # Net flux = downwelling - upwelling
    flux_sw_ref = rsd_ref - rsu_ref
    flux_sw     = rsd - rsu
    
    diff = np.max(np.abs(flux_sw - flux_sw_ref))
    print("Max. absolute error in SW net fluxes:")
    
    if (diff < 0.3):
         print("{:0.3f} = less than 0.3 W/m2 -------> Check PASSED!".format(diff))
    else:
         print("{:0.3f} = more than 0.3 W/m2 -------> Check NOT PASSED".format(diff))
    
    
    diff = mae(flux_sw, flux_sw_ref)
    print("Mean absolute error in SW net fluxes:")
    
    if (diff < 0.05):
         print("{:0.3f} = less than 0.05 W/m2 ------> Check PASSED!".format(diff))
    else:
         print("{:0.3f} = more than 0.05 W/m2 ------> Check NOT PASSED".format(diff))
    
    
    diff = mae(flux_sw[:,-1:],flux_sw_ref[:,-1:])
    print("Mean absolute error in downwelling SW flux at surface:")
    if (diff < 0.05):
         print("{:0.3f} = less than 0.05 W/m2 ------> Check PASSED!".format(diff))
    else:
         print("{:0.3f} = more than 0.05 W/m2 ------> Check NOT PASSED".format(diff))
    
    if do_sw_clear:
        flux_sw_clear_ref = rsd_clear_ref - rsu_clear_ref
        flux_sw_clear     = rsd_clear - rsu_clear
    
        diff = np.max(np.abs(flux_sw_clear-flux_sw_clear_ref))
    
        print("Max error in clear-sky shortwave net fluxes:")
        print("{:0.3f}".format(diff))
        
    print("----- HEATING RATES ------- ")
    
    diff = bias(dTdt_sw, dTdt_sw_ref)
    print("Bias in SW heating rate:")
    print("{:0.4f}".format(diff))

    inds_trop = pres_lay > 10000.0
    inds_strat = pres_lay < 10000.0

    diff = bias(dTdt_sw[inds_trop], dTdt_sw_ref[inds_trop])
    print("Bias in SW heating rate (troposphere):")
    print("{:0.4f}".format(diff))
    
    diff = mae(dTdt_sw, dTdt_sw_ref)
    
    print("Mean absolute error in SW heating rate:")
    
    if (diff < 0.02):
          print("{:0.3f} = less than 0.02 K/day ------> Check PASSED!".format(diff))
    else:
          print("{:0.3f} = more than 0.02 K/day ------> Check NOT PASSED".format(diff))
    
    diff = np.max(np.abs(dTdt_sw[inds_trop]-dTdt_sw_ref[inds_trop]))
    print("Max. absolute error in SW heating rates (troposphere):")
    if (diff < 0.2):
         print("{:0.3f} = less than 0.2 K/day ------> Check PASSED!".format(diff))
    else:
         print("{:0.3f} = more than 0.2 K/day ------> Check NOT PASSED".format(diff))



