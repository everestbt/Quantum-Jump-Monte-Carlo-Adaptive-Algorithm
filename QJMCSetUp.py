#Edited 12/3/17 Ben Everest
#Functions used to set up the algorithm and perform checks on the given
#variables
import scipy
import numpy
import sys
import random

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

	if (not isinstance(psi0,numpy.ndarray)):
		sys.exit("ERROR: the initial state is not a numpy ndarray")

#TESTED
def jumpOperatorsPaired(jumpOps):
	jumpOpsPaired = []
	for jumpOp in jumpOps:
		jumpOpsPaired.append(jumpOp.conjugate().transpose().dot(jumpOp))
	return jumpOpsPaired

#TESTED
def addExpectationSquared(eOps):
	for i in range(len(eOps)):
		eOps.append(eOps[i]*eOps[i])


def eResultsProduction(eOps,numberOfPoints):
	eResults = []
	for _ in range(len(eOps)):
		eResults.append(numpy.zeros(numberOfPoints))
	return eResults

def histogramProduction(histogramOptions,numberOfPoints):
	histograms = []
	for hist in histogramOptions:
		histograms.append(numpy.zeros((numberOfPoints,hist.numberOfBars)))

	return histograms

def HEffProduction(H, jumpOpsPaired):
	j=complex(0,1)
	HEff = H

	for jOpPaired in jumpOpsPaired:
		HEff = HEff - (j/2)*jOpPaired
	return HEff

def HEffExponentProduction(HEff, dt):
	j=complex(0,1)
	return scipy.linalg.expm(HEff.multiply(-j*dt))

#TODO make this use the HEffExponentProduction function
def HEffExponentSetProduction(H,jumpOpsPaired, dt, accuracyMagnitude):
	j=complex(0,1)

	HEff = HEffProduction(H, jumpOpsPaired)

	HEffExponentDt = scipy.linalg.expm(HEff.multiply(-j*dt))

	#Defines a HEff Exponent using a small time-steps
	HEffExponentDtSet = []
	for i in range(1,accuracyMagnitude + 1):
		HEffExponentDtSet.append(scipy.linalg.expm(
			HEff.multiply(-j*(dt/(pow(10,i))))))

	return HEffExponentDt, HEffExponentDtSet

#TESTED
def HEffExponentSetProductionBinary(H, jumpOpsPaired, deltaT, settings):
	HEff = HEffProduction(H, jumpOpsPaired)

	#Defines a HEff Exponent using a small time-steps
	HEffExponentDtSet = []
	dtSet = []
	dt = deltaT * 2
	#TODO add a safety check that the smallest dt in the list isn't smaller (if so just do the one)
	while (dt > settings.smallestDt):
		dt = dt/2
		dtSet.append(dt)
		HEffExponentDtSet.append(HEffExponentProduction(HEff, dt))

	return HEffExponentDtSet, dtSet

#TESTED
def randomInitialState(H):
	dim = H.get_shape()[0]
	psi0 = numpy.ndarray(shape=(dim,1),dtype=complex)
	for i in range(dim):
		r = random.random()
		u = random.random()
		psi0[i] = complex(r,u)

	psi0 = QJMCMath.normalise(psi0)
	return psi0
