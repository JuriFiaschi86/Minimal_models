#! /usr/bin/python3

# This script uses the Casas-Ibarra parametrization
# to transform the constrains and measurements on neutrino observables
# like mass splitting and mixing angles into constrains on the couplings of the model
# see arXiv:hep-ph/010365



import math
import numpy as np

########################
####### INPUT PARAMETERS
########################

### Input from SM
v = 2.46220569E+02

#### Input from specific BSM model (test point)
## Couplings:
#lambda4Input = 0.01
#lambda5Input = 0.055
### Masses:
## Fermion sector:
#mCPsiInput = 500
#mpsipsipInput = 920
## Scalar sector:
#MPHI2 = ((9.0000000E+04, 0), (0, 9.0000000E+04))

### Input from specific BSM model (realistic point)
# Couplings:
LAM1IN = ((0, 0), (0, 0))
lambda4Input = 5.0E-03
lambda5Input = 5.0E-03
## Masses:
# Fermion sector:
mCPsiInput = 1.0E+03
mpsipsipInput = 1.5E+03 
# Scalar sector:
MPHI2 = (((1.4E+03)**2, 0), (0, (3.0E+03)**2))
#MPHI2 = (((1.4E+03)**2, (2.5E+02)**2), ((2.5E+02)**2, (3.0E+03)**2))


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

######################################
####### LOOP GENERATED NEUTRINO MASSES
######################################

# Fermion mass matrix
M0 = np.matrix([[mCPsiInput, v*lambda5Input/math.sqrt(2), v*lambda4Input/math.sqrt(2)], [v*lambda5Input/math.sqrt(2), 0, mpsipsipInput], [v*lambda4Input/math.sqrt(2), mpsipsipInput, 0]])

# Eigenvalues and diagonalisation
eigenM0, UCHI = np.linalg.eig(M0)
#eigenM0 = abs(eigenM0)

# Scalar mass mass matrix

S0 = np.matrix([[MPHI2[0][0], MPHI2[0][1]], [MPHI2[1][0], MPHI2[1][1]]]) + (v**2) * np.matrix([[LAM1IN[0][0], LAM1IN[0][1]], [LAM1IN[1][0], LAM1IN[1][1]]])
eigenS0, US = np.linalg.eig(S0)
eigenS0 = abs(eigenS0)


## Test the diagonalisation fermions
##diagM0 = np.dot(UCHI, np.dot(M0, UCHI.transpose()))
#diagM0 = np.dot(UCHI.transpose(), np.dot(M0, UCHI))
#print(diagM0)
############### Note that Python gives the transpose matrix with respect to what one would expect. To fix the issue in the formula for M0nu I have to take the column of UCHI instead of the row

##Test the diagonalisation scalars
##diagS0 = np.dot(US, np.dot(S0, US.transpose()))
#diagS0 = np.dot(US.transpose(), np.dot(S0, US))
#print(diagS0)
############### Note that Python gives the transpose matrix with respect to what one would expect. To fix the issue in the formula for M0nu I have to take the column of UCHI instead of the row


# Model dependent loop integral for neutrino masses generation
############### Fix the issue of transpose of diagonalizing matrix taking the column of UCHI instead of the row

n_generation = 2
M0nu = [[None for j in range(n_generation)] for i in range(n_generation)]

Al = []
for i in range(0, n_generation):
    Al.append(sum((np.conjugate(np.array(UCHI[2,:]).flatten())**2)*(eigenM0**3 / (eigenS0[i] - eigenM0**2))*[math.log(x**2 / eigenS0[i]) for x in eigenM0]))

for i in range(0, n_generation):
    for j in range(0, n_generation):
        M0nu[i][j] = (1/(32*math.pi**2))*sum(Al*np.array(US[i,:]).flatten()*np.array(US[j,:]).flatten())

# Diagonalization of M0nu
eigenM0nu, UM = np.linalg.eig(M0nu)

DM_minusonehalf_item = []
for i in range(0, n_generation):
    DM_minusonehalf_item.append(1/math.sqrt(abs(eigenM0nu[i])))

DM_minusonehalf = np.diag(DM_minusonehalf_item)

#####################################
####### FUNCTION FOR PARAMETERS ARRAY
#####################################

def lambda6(theta):
    R = np.matrix([[0, math.cos(theta), -math.sin(theta)], [0, math.sin(theta), math.cos(theta)]])
    coupling = np.dot(UM, np.dot(DM_minusonehalf, (np.dot(R, np.dot(Dnusqrt, PMNS.H)))))
    return coupling.transpose()



print(lambda6(0))

###Test: I check if using these couplings I obtain the correct neutrino mass matrix (diagonalized)
#DM = np.diag(eigenM0nu)
#Dnu_test = np.dot(PMNS.H, np.dot(lambda6(0), np.dot(np.asmatrix(UM), np.dot(DM, np.dot(np.asmatrix(UM).H, np.dot(lambda6(0).transpose(), PMNS))))))

#print(Dnu_test)


#print(np.conjugate(np.array(UCHI[2,:]).flatten())**2)
#print((eigenM0**3 / (eigenS0[i] - eigenM0**2))*[math.log(x**2 / eigenS0[i]) for x in eigenM0])
#print(sum((np.conjugate(np.array(UCHI[2,:]).flatten())**2)*(eigenM0**3 / (eigenS0[i] - eigenM0**2))*[math.log(x**2 / eigenS0[i]) for x in eigenM0]))
