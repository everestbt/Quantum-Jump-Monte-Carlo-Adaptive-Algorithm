#Edited 12/3/17 Ben Everest
#An example of a 2-level system with plotting included
import sys
import qutip
import scipy
import numpy
import matplotlib.pyplot as plt

import QJMCAA
import QJMCMath
import QJMCMeasure

class Parameters:
	def __init__(self):
		self.omega = 1.0
		self.kappa = 1.0
		self.gamma = 1.0

#Define your initial state here using qutip commands
#This example sets a chain of spins to the up state
def initialStateDefine():
	#To define it using qutip
	psi0 = qutip.basis(2, 1)
	psi0 = psi0.full()

	#To define it manually
	#psi0 = numpy.ndarray(shape=(2,1),dtype=complex)
	#psi0[0] = complex(0,0)
	#psi0[1] = complex(1,0)
	#psi0 = QJMCMath.normalise(psi0) #just in case

	return psi0

#Defines the Hamiltonian
def Hamiltonian(parameters):
	sx = qutip.sigmax()

	H = parameters.omega * sx

	#All objects must be in a scipy sparse format
	H = H.full()

	#To define it manually
	#H = numpy.ndarray(shape=(2,2),dtype=complex)
	#H[0,1] = complex(parameters.omega,0)
	#H[1,0] = complex(parameters.omega,0)

	H = scipy.sparse.csc_matrix(H)
	return H

#Defines the jump operators
def jumpOperators(parameters):
	sp = qutip.sigmap()
	sm = qutip.sigmam()
	no = sp*sm

	jumpOps = []

	#Decay
	if parameters.kappa>0:
		jumpOps.append(numpy.sqrt(parameters.kappa) * sm)

	#Dephasing
	if parameters.gamma>0:
		jumpOps.append(numpy.sqrt(parameters.gamma) * no)

	for i in range(len(jumpOps)):
		jumpOps[i] = jumpOps[i].full()
		jumpOps[i]= scipy.sparse.csc_matrix(jumpOps[i])

	return jumpOps

#Defines the expectation operators that you want measured
def expectationOperators():
	si = qutip.qeye(2)
	sp = qutip.sigmap()
	sm = qutip.sigmam()
	no = sp*sm

	eOps = []

	eOps.append(no)


	for i in range(len(eOps)):
		eOps[i] = eOps[i].full()
		eOps[i] = scipy.sparse.csc_matrix(eOps[i])

	return eOps

def main(argv):
	#Sets up the parameters
	parameters = Parameters()
	parameters.omega = 1.0# you can also get it from terminal input float(argv[1])
	parameters.kappa = 1.0
	parameters.gamma = 0.0

	#Simulation settings
	settings = QJMCAA.Settings()
	settings.numberOfTrajectories = 100
	settings.T = 10.0
	settings.accuracyMagnitude = 5

	#Sets up how the data will be saved
	savingSettings = QJMCAA.SavingSettings()
	savingSettings.modelName = 'TwoLevel'
	#Adds the parameters to save by
	savingSettings.addParameter('Omega',parameters.omega)
	savingSettings.addParameter('Kappa',parameters.kappa)
	savingSettings.addParameter('Gamma',parameters.gamma)
	#Puts a name to the expectation values
	savingSettings.addExpectionSave('Population')
	#Sets the decimal places that we want the results to
	savingSettings.dataDecimalPlaces = 6
	#Turns on histograms
	savingSettings.histograms = True
	#You need to add each histogram you want, with the min, max and number of bars
	savingSettings.addHistogram(0.0,1.0,100)

	#Gets the defined Hamiltonian
	H = Hamiltonian(parameters)
	#Gets the defined jump operators
	jumpOps = jumpOperators(parameters)
	#Gets the defined expectation operators
	eOps = expectationOperators()

	#Gets the initial state you want to run it for
	psi0 = initialStateDefine()
	#This is the option to use a random state (have to pass a valid psi0 but it will be randomised)
	settings.randomInitialStates = False

	#Runs the simulation
	QJMCAA.QJMCRun(settings, savingSettings, H, jumpOps, eOps, psi0)

	#Plots the data
	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111)
	x = numpy.loadtxt(QJMCMeasure.nameTheFile(savingSettings,0))
	t = x[:,0]
	pop = x[:,1]
	var = x[:,2]
	ax1.plot(t,pop,'b')

	#This is how you add upper and lower boundaries (not very useful here)
	#upper = numpy.zeros(len(pop))
	#lower = numpy.zeros(len(pop))
	#for i in range(len(pop)):
		#upper[i] = pop[i] + numpy.sqrt(var[i])
		#lower[i] = pop[i] - numpy.sqrt(var[i])
	#ax1.plot(t,upper,'k--')
	#ax1.plot(t,lower,'k--')

	#Gets a solution from qutip master equation solution
	result = qutip.mesolve(parameters.omega *qutip.sigmax(),
		qutip.basis(2, 1),
		numpy.linspace(0.0, settings.T, settings.numberOfPoints),
		[numpy.sqrt(parameters.kappa) *qutip.sigmam(),
		numpy.sqrt(parameters.gamma) *qutip.sigmap()*qutip.sigmam()],
		[qutip.sigmap()*qutip.sigmam()])

	#Compares the two
	ax1.plot(result.times, result.expect[0],'r')
	ax1.set_xlabel('$t$')
	ax1.set_ylabel('$n$')
	ax1.autoscale(tight = 'True')

	#Plots the histogram data (note time is not stored in the histograms)
	name = QJMCMeasure.nameTheFile(savingSettings,0)
	name = name[:-4]
	name += 'Histogram'
	name += '.txt'
	hist = numpy.loadtxt(name)

	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111)

	popVals = numpy.arange(0,0.995,0.01)
	colorPlot = ax2.pcolor(popVals,t,hist,cmap = 'CMRmap')
	fig2.colorbar(colorPlot)

	ax2.set_xlabel('$n$')
	ax2.set_ylabel('$t$')
	ax2.autoscale(tight = 'True')

	plt.show()

if __name__ == "__main__":
	main(sys.argv[1:])
