#! /bin/python
#################################################################################################################################
#####################################################
# QMTIS QMTIMS QMTIS QMTIS QMTIS QMTIMS QMTIS QMTIS #
#####################################################
# Version 1.2
# Datae_ 2015-04-08
#####################################################
#
# Tool: Sys_Cstr (System Construction) ##############
#
# Developed by: Pengfei Ji 
#####################################################
# Note: This Python program is to be used to estimate the minimum thickness of each FMD zone, so that the Numerical stablity ####
# is guarantoeed. #################################################################################################################

import math

Delta_t_FDM = 0.05E-15

A = 3.57E06
B = 1.12E+11

T_e = 300.0
T_l = 300.0

V_Fermi = 1.39E+06

Delta_x_FDM = math.sqrt(Delta_t_FDM*(V_Fermi)**2/(A*T_e**2 + B*T_l))

print Delta_x_FDM
