# -*- coding: utf-8 -*-
from __future__ import division, print_function
from six.moves import range
import math
import numpy

# in GeV
HIGGS_VEV = 244.74
# For data: see PDG 2017
# assume NH
# Δm₂₁² = (7.53 ± 0.18) * 1e-5 eV²
# Δm₃₂² = (2.44 ± 0.06) * 1e-3 eV² (NH)
# Δm₃₂² = (2.51 ± 0.06) * 1e-3 eV² (IH)
# masses in GeV
M_NU_1 = 0
M_NU_2 = math.sqrt(7.53e-5) * 1e-9
M_NU_3 = math.sqrt(2.44e-21 + M_NU_2**2)
M_NU_2_LOW = math.sqrt(7.53e-5 - 0.18e-5) * 1e-9
M_NU_3_LOW = math.sqrt(2.44e-21 - 0.06e-21 + M_NU_2_LOW**2)
M_NU_2_HIGH = math.sqrt(7.53e-5 + 0.18e-5) * 1e-9
M_NU_3_HIGH = math.sqrt(2.44e-21 + 0.06e-21 + M_NU_2_LOW**2)
M_NU = numpy.array([M_NU_1, M_NU_2, M_NU_3])
D_NU_onehalf = numpy.matrix(numpy.diag(M_NU**(1/2)))
#D_NU_minusonehalf = numpy.matrix(numpy.diag([0] + list(M_NU[1:]**(1/2))))

# PMNS matrix
# sin(θ₁₂)² = 0.304 ± 0.014
# sin(θ₁₃)² = (2.19 ± 0.12) * 1e-2
# sin(θ₂₃)² = 0.51 ± 0.05 (NH)
# sin(θ₂₃)² = 0.50 ± 0.05 (IH)
THETA_12 = math.asin(math.sqrt(0.304))
THETA_13 = math.asin(math.sqrt(2.19e-2))
THETA_23 = math.asin(math.sqrt(0.51))
C12 = math.cos(THETA_12)
C13 = math.cos(THETA_13)
C23 = math.cos(THETA_23)
S12 = math.sin(THETA_12)
S13 = math.sin(THETA_13)
S23 = math.sin(THETA_23)
PMNS = numpy.matrix([
	[C12 * C13, S12 * C13, S13],
	[-S12 * C23 - C12 * S23 * S13, C12 * C23 - S12 * S23 * S13, S23 * C13],
	[S12 * S23 - C12 * C23 * S13, -C12 * S23 - S12 * C23 * S13, C23 * C13],
])

def T13Balpha02g(model, theta):
	N_FERMION = 3
	N_SCALAR = 2
	v = HIGGS_VEV
	# parameter values
	l4 = model.param_value_all('lambda4')
	l5 = model.param_value_all('lambda5')
	mP = model.param_value_all('mCPsi')
	mpp = model.param_value_all('mpsipsip')
	l1 = numpy.matrix(model.param_value_all('lambda1'))
	mphi2 = numpy.matrix(model.param_value_all('mphi2'))
	# fermion and scalar mass matrices
	M_chi = numpy.matrix([
		[mP, v * l5/math.sqrt(2), v * l4/math.sqrt(2)],
		[v * l5/math.sqrt(2), 0, mpp],
		[v * l4/math.sqrt(2), mpp, 0],
	])
	M_eta2 = mphi2 + v**2 * l1
	# diagonalize mass matrices
	m_chi, U_chi_T = numpy.linalg.eig(M_chi)
	m_eta2, O_eta_T = numpy.linalg.eig(M_eta2)
	U_chi = numpy.matrix(U_chi_T.transpose())
	O_eta = numpy.matrix(O_eta_T.transpose())
	# neutrino mass matrix
	M = numpy.empty((N_SCALAR, N_SCALAR))
	for m in range(N_SCALAR):
		for n in range(N_SCALAR):
			M[m, n] = 1/(32 * math.pi**2) * sum(
				O_eta[l, m] * O_eta[l, n]
					* U_chi[k, N_FERMION - 1].conjugate()**2
					* m_chi[k]**3/(m_eta2[l] - m_chi[k]**2)
					* math.log(m_chi[k]**2/m_eta2[l])
				for k in range(N_FERMION)
				for l in range(N_SCALAR)
			)
	D_M, U_M_T = numpy.linalg.eig(M)
	print('M')
	print(M)
	print(D_M)
	assert all(D_M != 0), 'D_M singular!'
	D_M = numpy.abs(D_M)
#	D_M_onehalf = numpy.matrix(numpy.diag(D_M**(1/2)))
	D_M_minusonehalf = numpy.matrix(numpy.diag(D_M**(-1/2)))
	U_M = numpy.matrix(U_M_T.transpose())
	# matrix R
	R = numpy.matrix([
		[0, 0],
		[math.cos(theta), math.sin(theta)],
		[-math.sin(theta), math.cos(theta)],
	])
	# result
	l6 = (
		PMNS.conjugate()
			* D_NU_onehalf * R * D_M_minusonehalf
		* U_M.conjugate()
	)
#	l6 = numpy.matrix(model.param_value_all('lambda6'))
	print('l6')
	print(PMNS)
	print(D_NU_onehalf)
	print(R)
	print(D_M_minusonehalf)
	print(U_M)
	print(l6)
	# M_nu
	print('M_nu')
	M_nu = numpy.empty((3, 3))
	for i in range(3):
		for j in range(3):
			M_nu[i, j] = sum(l6[i, m] * l6[j, n] * M[m, n]
				for m in range(N_SCALAR)
				for n in range(N_SCALAR)
			)
	print(M_nu)
	D_M_nu, U_nu_T = numpy.linalg.eig(M_nu)
	print(D_M_nu)
	return 'lambda6', l6

