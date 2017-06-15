#Edited 13/3/17 Ben Everest
#The functions used to perform a jump in the low memory trajectory
import random

import QJMCMath
import QJMCJump
import QJMCSetUp

def approxFindDt(tStart, deltaT, smallestDt, HEff, psi, r):
	t = tStart + deltaT
	previousPsi = psi
	dt = deltaT
	while dt > smallestDt:
		psi = previousPsi
		t -= dt
		survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)
		dt = dt / 10
		#TODO this is not used but need to check why
		#HEffExponentDt = QJMCSetUp.HEffExponentProduction(HEff, dt)
		while survivialProbability>r:
			t += dt
			survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)

	return t, psi

#Catches up to the correct time
def catchUpApprox(t, tEnd, smallestDt, psi, jumpOps, jumpOpsPaired,HEff):
	#Generates the random number for the new time selection
	r = random.random()
	#TODO change this to a do a close check
	while (t < tEnd):
		#Evolves to the appropriate time
		dt = tEnd - t
		HEffExponentDt = QJMCSetUp.HEffExponentProduction(HEff, dt)
		previousPsi = psi
		psi = HEffExponentDt.dot(psi)
		#Progress the time
		t += dt
		#Checks if it has jumped
		survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)
		if survivialProbability<r:
			tStart = t - dt
			t, psi = jump(tStart, dt, smallestDt, previousPsi, r, jumpOps,
				jumpOpsPaired, HEff)
			r = random.random()

	return t, psi, r

def jump(tStart, deltaT, smallestDt, psi, r, jumpOps, jumpOpsPaired,
		HEff):
	#Extracts where the jump occured
	if deltaT > smallestDt:
		t, psi = approxFindDt(tStart, deltaT, smallestDt, HEff, psi, r)
	else:
		HEffExponentDt = HEffExponentProduction(HEff, deltaT)
		psi = HEffExponentDt.dot(psi)
		t = tStart + deltaT

	#Selects which jump occurs
	jumpSelect = QJMCJump.jumpSelection(psi, jumpOpsPaired)
	#Performs the jump
	psi = jumpOps[jumpSelect].dot(psi)
	psi = QJMCMath.normalise(psi)

	return t, psi
