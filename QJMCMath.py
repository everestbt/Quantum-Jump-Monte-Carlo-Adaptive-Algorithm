import numpy
import math

#TESTED
def calculateSquareOfWavefunction(psi):
	P = numpy.asscalar(numpy.real(numpy.dot(numpy.conjugate(numpy.transpose(psi)),psi)))
	return P

#TESTED
def neighbour(pos, N, shift):
	return (pos + N + shift)%N

#TESTED
def normalise(psi):
	P = calculateSquareOfWavefunction(psi)
	psiN = psi / numpy.sqrt(P)
	return psiN

def magnitude(x):
	return int(math.log10(x))