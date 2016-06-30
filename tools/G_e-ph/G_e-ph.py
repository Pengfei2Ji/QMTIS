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
# Note: This Python program is to be used to caluclate the Electron Temperature depdendent Electron-Phonon Coupling Factor #######
from math import exp
import numpy as np

Conv_Factor_Ha2eV = 27.2114								# 1 Ha = 27.2114 eV = 2625.4996 kJ/mol
Conv_Factor_eV2J = 1.6*10**(-19)/(((4.2086734254742355*10**(-10))**3)/4.0)		# From eV/Volume_Atom to J/m^3

Pi = 3.1415
h = 6.62607E-34
h_Bar = 6.582119514E-16									# Unit: eV*s
kB = 8.6173324*10**(-5)         							# Unit: eV/K

Lamda_Omega_Square = 0.0		
Lamda_Omega_Square_eV = 0.0
Lamda = 0.0

###################################################################################################################################
# Update Field === Update Field #
g_E_Fermi = 7.45843671/Conv_Factor_Ha2eV						# Fermi energy is: 0.28121395
Mu = -0.14894321*Conv_Factor_Ha2eV							# Unit: eV
Te = 25000.00										# Unit: K
dOmega = 2.2681391872275680E-006
####################################################################################################################################
G_Electron_Phonon = 0.0									# Unit: W/(m^3*K)
Calclus = 0.0
dE = 0.00050*Conv_Factor_Ha2eV

DOS_File = open ('Ag_PAW_Laser_out_DOS', 'r')

for line in DOS_File:
	Col_DOS = line.split()
	E = float(Col_DOS[0])*Conv_Factor_Ha2eV
	g = float(Col_DOS[1])/Conv_Factor_Ha2eV
	f = 1/(np.exp((E-Mu)/(kB*Te)) + 1)
#	df = ((np.exp((E-Mu)/(kB*Te)) + 1)**(-2))*np.exp((E-Mu)/(kB*Te))*(-Mu)/(kB*Te)
        f1 = 1/(np.exp(((E-dE)-Mu)/(kB*Te)) + 1)
        f2 = 1/(np.exp(((E+dE)-Mu)/(kB*Te)) + 1)
	df = (f2-f1)/(2*dE)
	Calclus += -g**2*df*dE
print "The Calclus result is: %f." % Calclus
####################################################################################################################################
A2F_Dat = open ('teph_5.out_ep_A2F', 'r')

for line in A2F_Dat:
	Col_A2F_Dat = line.split()
	Lamda_Omega_Square_eV += 2*float(Col_A2F_Dat[0])*float(Col_A2F_Dat[1])*((27.2114)**2)*10**6*dOmega
	Lamda_Omega_Square += 2*float(Col_A2F_Dat[0])*float(Col_A2F_Dat[1])*((6.57969e15)**2)*dOmega
	Lamda += 2*float(Col_A2F_Dat[1])/float(Col_A2F_Dat[0])*dOmega
print "Lamda is: %f" % Lamda
print "The calculated Lamda_Omega_Square is: %5e" % Lamda_Omega_Square
print "The calculated Lamda_Omega_Square in mev^2 is: %f" % Lamda_Omega_Square_eV
####################################################################################################################################
G_Electron_Phonon = Pi*h_Bar*kB*Lamda_Omega_Square/g_E_Fermi*Calclus*Conv_Factor_Ha2eV*Conv_Factor_eV2J
print "The electron-phonon coupling factor is: %5e" % G_Electron_Phonon
