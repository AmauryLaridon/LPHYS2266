#######################################################################################################################################
############################################### Personnal Exercice Computation ########################################################
# Author: Amaury Laridon 
# Course: LPHYS2266 - Physics of the Upper Atmosphere and Space 
# Goal: Computation of the Pannekoek-Rosseland electric field Intensity and the ionic density scale height 
#       at the exobase assuming protonosphere.
# Date: 14/04/23
#######################################################################################################################################
#######################################################################################################################################
##### Physical constant #####

m_p = 1.672649 * (10^(-27))
m_e = 9.109 * (10^(-31))
e = -1.602 * (10^(-19))
k = 1.38 * (10^(-23))

##### Parameters #####

g_0_ground = 9.81
r_0 = 6371000
h_exobase = 2000000
r_c = r_0 + h_exobase
g_0_exobase = g_0_ground * ((r_0 / r_c)^2)
T_0_ground = 288.15
T_0_exobase = 1000

##################################### Computation of the Pannekoek-Rosseland Electric Field Intensity ####################################
num1 = (m_p - m_e)
denum1 = 2 * (e)
num2 = r_0
denum2 = r_c

E = g_0_ground * (num1 / denum1) * ((num2 / denum2)^2)
println("Intensity of the Pannekoek-Rosseland Electric Field at the exobase = ", E, " N/C")

################################################ Computation of the ionic density scale height ############################################

num3 = k * T_0_exobase
denum3 = (m_p) * (g_0_exobase)
H = 2 * (num3 / denum3) # [m]

println(g_0_exobase)
println("Ionic density scale height = ", H / 1000, "km")
