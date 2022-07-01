#! /usr/bin/python3

import math
import numpy as np

########################
####### INPUT PARAMETERS
########################

### Input from SM
v = 2.46220569E+02


### Input from experimental constrains
## Neutrino mixing:
th12 = 33.62 * math.pi / 180
th23 = 47.2 * math.pi / 180
th13 = 8.54 * math.pi / 180
#deltaCP = 234 * math.pi / 180
deltaCP = 0

PMNS23 = np.matrix([[1, 0, 0], [0, math.cos(th23), math.sin(th23)], [0, -math.sin(th23), math.cos(th23)]])
#PMNS13 = np.matrix([[math.cos(th13), 0, math.sin(th13)*cmath.exp(-deltaCP j)], [0, 1, 0], [-math.sin(th13)*cmath.exp(+deltaCP j), 0, math.cos(th13)]])
PMNS13 = np.matrix([[math.cos(th13), 0, math.sin(th13)], [0, 1, 0], [-math.sin(th13), 0, math.cos(th13)]])
PMNS12 = np.matrix([[math.cos(th12), math.sin(th12), 0], [-math.sin(th12), math.cos(th12), 0], [0, 0, 1]])

PMNS = np.dot(PMNS23, np.dot(PMNS13, PMNS12))

## Neutrino mass splitting
mnu1 = 0
mnu2 = math.sqrt(7.53E-05)*1.0E-09
mnu3 = math.sqrt(2.44E-03 + mnu2**2)*1.0E-09

Dnu = np.matrix([[mnu1, 0, 0], [0, mnu2, 0], [0, 0, mnu3]])
Dnusqrt = np.matrix([[math.sqrt(mnu1), 0, 0], [0, math.sqrt(mnu2), 0], [0, 0, math.sqrt(mnu3)]])


####################
### MODEL PARAMETERS
####################

M = [1.0E+00, 1.0E+06, 1.0E+06]

eps = 1.0E-04
m4 = math.sqrt(1.0E-05)


lambda4Input = 5
lambda5Input = 5
lambda6Input = 5


######################
### DERIVED QUANTITIES
######################

m2 = math.sqrt((lambda4Input + lambda5Input)*(v**2)/2)
m3 = math.sqrt(lambda6Input*(v**2)/2)

A = m2**2 + (lambda4Input + lambda5Input)*(v**2)/2
CR = m3**2 + lambda6Input*(v**2)/2 + m4**2
CI = m3**2 + lambda6Input*(v**2)/2 - m4**2

mu = math.sqrt(2*(A*CR - eps*(A + CR)) / (v**2))

B = mu*v/math.sqrt(2)


#################
### SCALAR MASSES
#################


MR2 = np.matrix([[A, B], [B, CR]])
MI2 = np.matrix([[A, B], [B, CI]])

#MR = np.matrix([[math.sqrt(A), math.sqrt(B)], [math.sqrt(B), math.sqrt(CR)]])
#MI = np.matrix([[math.sqrt(A), math.sqrt(B)], [math.sqrt(B), math.sqrt(CI)]])


eigenMR, UR = np.linalg.eig(MR2)
eigenMI, UI = np.linalg.eig(MI2)


##################################
### LOOP GENERATED NEUTRINO MASSES
##################################

#M0nu = (1/(16*(math.pi**2)))*sum(m*(((UR[0,0]**2)*(eigenMR[0]**2) / (eigenMR[0]**2 - m**2) * math.log((eigenMR[0]**2)/(m**2))) + ((UR[1,0]**2)*(eigenMR[1]**2) / (eigenMR[1]**2 - m**2) * math.log((eigenMR[1]**2)/(m**2))) - ((UI[0,0]**2)*(eigenMI[0]**2) / (eigenMI[0]**2 - m**2) * math.log((eigenMI[0]**2)/(m**2))) - ((UI[1,0]**2)*(eigenMI[1]**2) / (eigenMI[1]**2 - m**2) * math.log((eigenMI[1]**2)/(m**2)))) for m in M)


M0nu = abs((1/(16*(math.pi**2)))*sum(m*(((UR[0,0]**2)*(eigenMR[0]) / (eigenMR[0] - m**2) * math.log((eigenMR[0])/(m**2))) + ((UR[1,0]**2)*(eigenMR[1]) / (eigenMR[1] - m**2) * math.log((eigenMR[1])/(m**2))) - ((UI[0,0]**2)*(eigenMI[0]) / (eigenMI[0] - m**2) * math.log((eigenMI[0])/(m**2))) - ((UI[1,0]**2)*(eigenMI[1]) / (eigenMI[1] - m**2) * math.log((eigenMI[1])/(m**2)))) for m in M))


def lambda6(theta):
    R = np.matrix([[0, math.cos(theta), -math.sin(theta)], [0, math.sin(theta), math.cos(theta)]])
    coupling = pow(M0nu,-0.5)*(np.dot(R, np.dot(Dnusqrt, PMNS.H)))
    return coupling



print(lambda6(0))
