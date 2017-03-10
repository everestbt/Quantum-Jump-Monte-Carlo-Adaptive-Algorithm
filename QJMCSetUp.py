import scipy
import numpy
import sys
import random
import math

import QJMCAA
import QJMCMath

#TESTED
def dimensionTest(H,jumpOps,eOps,psi0):
	#Compare all to the hamiltonian
	dim = H.get_shape()

	c = 0
	for item in jumpOps:
		c+=1
		if (item.get_shape() != dim):
			sys.exit("ERROR: the "+ str(c) +" jump operator (or more) are the wrong dimension with respect to H")

	c = 0
	for item in eOps:
		c+=1
		if (item.get_shape() != dim):
			sys.exit("ERROR: the "+ str(c) +" jump operator (or more) are the wrong dimension with respect to H")

	if (psi0.shape[0] != dim[0]):
		sys.exit("ERROR: the initial state is the wrong dimension")

	if (psi0.shape[1] != 1):
		sys.exit("ERROR: the initial state is the wrong dimension")

#TESTED
def typeTest(settings, savingSettings, H, jumpOps, eOps, psi0):
	#Checks the settings data typeTest
	if (not isinstance(settings, QJMCAA.Settings)):
		sys.exit("ERROR: the settings object is not the correct class")

	if (not isinstance(savingSettings, QJMCAA.SavingSettings)):
		sys.exit("ERROR: the savingSettings object is not the correct class")

	if (not scipy.sparse.issparse(H)):
		sys.exit("ERROR: the H is not a sparse scipy array")

	c = 0
	for item in jumpOps:
		c+=1
		if (not scipy.sparse.issparse(item)):
			sys.exit("ERROR: the "+ str(c) +" jump operator (or more) is not a sparse scipy array")

	c = 0
	for item in eOps:
		c+=1
		if (not scipy.sparse.issparse(item)):
			sys.exit("ERROR: the "+ str(c) +" expectation operator (or more) is not a sparse scipy array")

	if (not type(psi0) is numpy.ndarray):
		sys.exit("ERROR: the initial state is not a numpy ndarray")

#TESTED
def jumpOperatorsPaired(jumpOps):
	jumpOpsPaired = []
	for i in range(len(jumpOps)):
		jumpOpsPaired.append(jumpOps[i].conjugate().transpose().dot(jumpOps[i]))
	return jumpOpsPaired


def addExpectationSquared(eOps):
	for i in range(len(eOps)):
		eOps.append(eOps[i]*eOps[i])


def eResultsProduction(eOps,numberOfPoints):
	eResults = []
	for i in range(len(eOps)):
		eResults.append(numpy.zeros(numberOfPoints))
	return eResults

def histogramProduction(histogramOptions,numberOfPoints):
	histograms = []
	for hist in histogramOptions:
		histograms.append(numpy.zeros((numberOfPoints,hist.numberOfBars)))

	return histograms

def HEffExponentProduction(H,jumpOpsPaired, dt, accuracyMagnitude):
	j=complex(0,1)

	HEff = H

	#This is now HEffSparse but we use just H to avoid doubling the memory
	for i in range(len(jumpOpsPaired)):
		HEff = HEff - (j/2)*jumpOpsPaired[i]

	HEffExponentDt = scipy.linalg.expm(HEff.multiply(-j*dt))

	#Defines a HEff Exponent using a small time-steps
	HEffExponentDtSet = []
	for i in range(1,accuracyMagnitude + 1):
		HEffExponentDtSet.append(scipy.linalg.expm(
			HEff.multiply(-j*(dt/(pow(10,i))))))

	return HEffExponentDt, HEffExponentDtSet

#TESTED
def randomInitialState(H):
	j = complex(0,1)

	dim = H.get_shape()[0]
	psi0 = numpy.ndarray(shape=(dim,1),dtype=complex)
	for i in range(dim):
		r = random.random()
		u = random.random()
		psi0[i] = complex(r,u)

	psi0 = QJMCMath.normalise(psi0)
	return psi0
