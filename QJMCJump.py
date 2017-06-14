#Edited 12/3/17 Ben Everest
#The functions used to perform a jump in the trajectory
import numpy
import random

import QJMCMath
import QJMCEvolve

def approxFindDt(tStart, deltaT, HEffExponentDtSet, psi, r):
	t = tStart
	previousPsi = psi
	dt = 0
	for i in range(1,len(HEffExponentDtSet)+1):
		psi = previousPsi
		t -= dt
		survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)
		dt = deltaT/(pow(10.0,i))
		while survivialProbability>r:
			previousPsi = psi
			psi = HEffExponentDtSet[i-1].dot(psi)
			t += dt
			survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)

	return t, psi

def jumpSelection(psi, jumpOpsPaired):
	if (len(jumpOpsPaired) == 1):
		return 0
	#Gets a new random number for jump selection
	r = random.random()
	#print(r)
	psiLeft = numpy.conjugate(numpy.transpose(psi))[0]
	P = numpy.zeros(len(jumpOpsPaired))
	SumPj = 0
	for i, pairedJumpOp in enumerate(jumpOpsPaired):
		P[i] = numpy.real(numpy.dot(psiLeft,pairedJumpOp.dot(psi))[0])
		SumPj += P[i]
	r *= SumPj
	jumpSelect = -1
	while r>0:
		jumpSelect += 1
		r -= P[jumpSelect]
	return jumpSelect

#Catches up to the correct time
def catchUpApprox(t, tEnd, deltaT, psi, jumpOps, jumpOpsPaired,HEffExponentDtSet):
	#Generates the random number for the new time selection
	r = random.random()
	for i in range(1,len(HEffExponentDtSet)+1):
		#Start the change in time at the largest first and then get smaller
		dt = deltaT/(pow(10.0,i))
		while ((tEnd - t) > dt):
			#Save the previous wave function
			prev = psi
			#Evolve it
			psi = HEffExponentDtSet[i-1].dot(psi)
			#Progress the time
			t += dt
			#Check if it has jumped
			survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)
			if survivialProbability<r:
				#Peform the jump in the last time step
				tStart = t - dt
				t, psi = jump(tStart, deltaT, prev, r,
					jumpOps, jumpOpsPaired, HEffExponentDtSet)
				#Produce a new random number for time selection
				r = random.random()
	return t, psi, r



def jump(tStart, deltaT, psi, r, jumpOps, jumpOpsPaired,
		HEffExponentDtSet):
	#Extracts where the jump occured
	t, psi = approxFindDt(tStart, deltaT,  HEffExponentDtSet, psi, r)
	#Selects which jump occurs
	jumpSelect = jumpSelection(psi, jumpOpsPaired)
	#Performs the jump
	psi = jumpOps[jumpSelect].dot(psi)
	psi = QJMCMath.normalise(psi)

	return t, psi

#Catches up to the correct time
#TODO finish this function
def catchUpApproxBinary(t, psi, dtSet, HEffExponentDtSet, catchUpSet):
	#Generates the random number for the new time selection
	r = random.random()
	catchUpReturn = catchUpSet
	for i in list(reversed(catchUpSet)):
		previousPsi = psi
		t, psi = QJMCEvolve.evolvePsiByExponent(psi, t, dtSet[i], HEffExponentDtSet[i])
		catchUpReturn.remove(i)
		survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)
		if survivialProbability<r:
			return t - dtSet[i], previousPsi, r, i+1, catchUpReturn
	return t, psi, r, 0, catchUpReturn

def approxFindDtBinary(t, HEffExponentDtSet, dtSet, psi, r, startingIndex, catchUpSet):
	previousPsi = psi
	#Does a binary search for the jump position
	for i in range(startingIndex,len(HEffExponentDtSet)-1):
		t, psi = QJMCEvolve.evolvePsiByExponent(psi, t, dtSet[i], HEffExponentDtSet[i])
		survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)
		if survivialProbability>r:
			previousPsi = psi
		else:
			t -= dtSet[i]
			psi = previousPsi
			catchUpSet.append(i)
	#Does the final search which needs to correct for a +ve result
	t, psi = QJMCEvolve.evolvePsiByExponent(psi, t, dtSet[-1], HEffExponentDtSet[-1])
	survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)
	if survivialProbability>r:
		t, psi = QJMCEvolve.evolvePsiByExponent(psi, t, dtSet[-1], HEffExponentDtSet[-1])
	else:
		catchUpSet.append(len(HEffExponentDtSet)-1)

	return t, psi, catchUpSet

#TESTED
def jumpBinary(tStart, psi, r, jumpOps, jumpOpsPaired, dtSet, HEffExponentDtSet):
	t = tStart
	startingIndex = 1
	catchUpSet = []
	while (tStart + dtSet[0] - t) > dtSet[-1]/2:
		#Extracts where the jump occured
		if (startingIndex != len(dtSet)):
			t, psi, catchUpSet = approxFindDtBinary(t, HEffExponentDtSet,
				dtSet, psi, r, startingIndex, catchUpSet)
		else:
			t, psi = QJMCEvolve.evolvePsiByExponent(psi, t, dtSet[-1],
				HEffExponentDtSet[-1])
		#print(t)
		#print(catchUpSet)
		#Selects which jump occurs
		jumpSelect = jumpSelection(psi, jumpOpsPaired)
		#Performs the jump
		psi = jumpOps[jumpSelect].dot(psi)
		psi = QJMCMath.normalise(psi)
		#Performs catch-up

		t, psi, r, startingIndex, catchUpSet = catchUpApproxBinary(t, psi,
			dtSet, HEffExponentDtSet, catchUpSet)

	return t, psi, r
