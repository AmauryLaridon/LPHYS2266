#######################################################################################################################################
################################################### Wind Mission Data Analysis ########################################################
# Author: Amaury Laridon
# Course: LPHYS2266 - Physics of the Upper Atmosphere and Space
# Goal: Extraction of some data from the Wind NASA mission and analysis.
#       The goal is to see wether the observed time variations in several variables are due to change
#       in the orbit of the satellite or to solar activity.
# Data: Wind data from NASA available at https://cdaweb.gsfc.nasa.gov/index.html
# Date: 12/04/23
#######################################################################################################################################
#######################################################################################################################################

######################################################## Librairies Loading ###########################################################
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import csv
######################################################## Data Extraction ###############################################################

### Extraction of the data from the CSV files to a DataFrame object ###

data_var = []
with open('/home/amaury/Bureau/LPHYS2266 - Physics of the upper atmosphere and space/Projet/Data Analysis/Data/OMNI_HRO_1MIN_126357.csv', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in my_reader:
        data_var.append(row)
print(np.shape(data_var))
time_var_str = data_var[:]
print(time_var_str)
print(np.shape(time_var_str))


""" data_var = DataFrame(CSV.File(
    "/home/amaury/Bureau/LPHYS2266 - Physics of the upper atmosphere and space/Projet/Data Analysis/Data/OMNI_HRO_1MIN_126357.csv"))
data_ind = DataFrame(CSV.File(
    "/home/amaury/Bureau/LPHYS2266 - Physics of the upper atmosphere and space/Projet/Data Analysis/Data/OMNI2_H0_MRG1HR_126357.csv"))
data_orb = DataFrame(CSV.File(
    "/home/amaury/Bureau/LPHYS2266 - Physics of the upper atmosphere and space/Projet/Data Analysis/Data/WI_OR_PRE_126357.csv"))
### Definition of the arrays ###

time_var_str = data_var[:, 1]
time_var = range(1, length(time_var_str))
flow_speed = data_var[:, 2]
prot_density = data_var[:, 3]
temperature = data_var[:, 4]
# println(length(flow_speed))

time_ind_str = data_ind[:, 1]
time_ind = range(1, length(time_ind_str))
Kp = data_ind[:, 2]
Dst = data_ind[:, 3]
Lyman_α = data_ind[:, 4]
# println(length(Kp))

time_orb_str = data_orb[:, 1]
time_orb = range(1, length(time_orb_str))
GSE_Z = data_orb[:, 4]
# println(length(GSE_Z)) """
"""
##################################################### Extrapolation of the data ########################################################

# We need to extrapolate some of the data since the series doesn't have the same length due to differences
# in temporal resolation of the instruments.

new_grid = time_var

interp_func_Kp = sc_interp.interp1d(time_ind, Kp, fill_value="extrapolate")
interp_Kp = interp_func_Kp(new_grid)

interp_func_Dst = sc_interp.interp1d(time_ind, Dst, fill_value="extrapolate")
interp_Dst = interp_func_Dst(new_grid)

interp_func_Lyman = sc_interp.interp1d(
    time_ind, Lyman_α, fill_value="extrapolate")
interp_Lyman = interp_func_Lyman(new_grid)

interp_func_GSE = sc_interp.interp1d(time_orb, GSE_Z, fill_value="extrapolate")
interp_GSE_Z = interp_func_GSE(new_grid)

println(length(interp_Dst))
println(length(interp_Lyman))
println(length(interp_GSE_Z))
############################################################ Display ##################################################################

data_plot = scatter(interp_Kp, flow_speed)
title!("Plasma Flow Speed - Kp index relation")
# plot!(new_grid, interp_Kp, label="Interp Kp")
# plot!(time_mod, regr_lin_opt, label ="Wieghted least-squares fit ")
# display(data_plot)
savefig("flow_speed(Kp).png")
"""
