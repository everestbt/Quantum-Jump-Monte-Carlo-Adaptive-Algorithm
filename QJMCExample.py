#Edited 12/3/17 Ben Everest
#An example of a lattice system
import sys
import qutip
import QJMCMath
import scipy
import numpy as np

import QJMCAA

class Parameters:
	def __init__(self):
		self.omega = 1.0
		self.kappa = 100.0
		self.gamma = 1.0

class Lattice:
	def __init__(self):
		self.numberOfSites = 1

#Define your initial state here using qutip commands
#This example sets a chain of spins to the up state
def initialStateDefine(N):
	psi_list = []
	for _ in range(N):
		psi_list.append(qutip.basis(2, 0))
	psi0 = qutip.tensor(psi_list)

	psi0 = psi0.full()
	return psi0

#Defines the Hamiltonian
def Hamiltonian(parameters, lattice):
	#Constructs the operators
	#(these are the essential operators that are needed for spin)
	si = qutip.qeye(2)
	sx = qutip.sigmax()
	#sy = qutip.sigmay()
	#sz = qutip.sigmaz()
	sp = qutip.sigmap()
	sm = qutip.sigmam()
	no = sp*sm

	#Constructs operators for the chain of length N
	#si_list = []
	sx_list = []
	#sy_list = []
	#sz_list = []
	#sp_list = []
	#sm_list = []
	no_list = []

	#Runs over each site defining the operator on that site
	for k in range(lattice.numberOfSites):
		#Puts an indentity on every site
		op_list = [si] * lattice.numberOfSites
		#Defines the sigma_x on site k
		op_list[k] = sx
		sx_list.append(qutip.tensor(op_list))
		#op_list[k] = sy
		#sy_list.append(qutip.tensor(op_list))
		#op_list[k] = sz
		#sz_list.append(qutip.tensor(op_list))
		#op_list[k] = sp
		#sp_list.append(qutip.tensor(op_list))
		#op_list[k] = sm
		#sm_list.append(qutip.tensor(op_list))
		op_list[k] = no
		no_list.append(qutip.tensor(op_list))

	#Constructs the Hamiltonian
	H = 0

	#Periodic boundary conditions
	#TODO define the periodic boundary conditions and closed in QJMCMath.neighbour operator with a choice
	for k in range(lattice.numberOfSites):
		H += (parameters.omega* (no_list[QJMCMath.neighbour(k, lattice.numberOfSites, -1)] +
			no_list[QJMCMath.neighbour(k, lattice.numberOfSites, 1)])*sx_list[k])

	#All objects must be in a scipy sparse format
	H = H.full()
	H = scipy.sparse.csc_matrix(H)
	return H

#Defines the jump operators
def jumpOperators(parameters, lattice):
	#Constructs the operators
	si = qutip.qeye(2)
	#sx = qutip.sigmax()
	#sy = qutip.sigmay()
	#sz = qutip.sigmaz()
	sp = qutip.sigmap()
	sm = qutip.sigmam()
	no = sp*sm
	#print no

	#Constructs operators for the chain of length N
	#si_list = []
	#sx_list = []
	#sy_list = []
	#sz_list = []
	sp_list = []
	sm_list = []
	no_list = []

	#si_list.append(qutip.tensor([si] * N))

	for k in range(lattice.numberOfSites):
		op_list = [si] * lattice.numberOfSites
		#op_list[k] = sx
		#sx_list.append(qutip.tensor(op_list))
		#op_list[k] = sy
		#sy_list.append(qutip.tensor(op_list))
		#op_list[k] = sz
		#sz_list.append(qutip.tensor(op_list))
		op_list[k] = sp
		sp_list.append(qutip.tensor(op_list))
		op_list[k] = sm
		sm_list.append(qutip.tensor(op_list))
		op_list[k] = no
		no_list.append(qutip.tensor(op_list))

	#Collapse operators
	c_op_list = []

	#Flips
	if parameters.kappa>0:
		#Flip down check right
		for k in range(lattice.numberOfSites):
			c_op_list.append(np.sqrt(parameters.kappa) * sm_list[k] *
				no_list[QJMCMath.neighbour(k, lattice.numberOfSites, 1)])
		#Flip down check left
		for k in range(lattice.numberOfSites):
			c_op_list.append(np.sqrt(parameters.kappa) * sm_list[k] *
				no_list[QJMCMath.neighbour(k, lattice.numberOfSites, -1)])
		#Flip up check right
		for k in range(lattice.numberOfSites):
			c_op_list.append(np.sqrt(parameters.kappa) * sp_list[k] *
				no_list[QJMCMath.neighbour(k, lattice.numberOfSites, 1)])
		#Flip up check left
		for k in range(lattice.numberOfSites):
			c_op_list.append(np.sqrt(parameters.kappa) * sp_list[k] *
				no_list[QJMCMath.neighbour(k, lattice.numberOfSites, -1)])

	#Decay (the greater than 0 prevents the production of empty jump opertors)
	if parameters.gamma >0:
		for k in range(lattice.numberOfSites):
			c_op_list.append(np.sqrt(parameters.gamma) * sm_list[k])

	#Converts the existing arrays into sparse arrays
	for i in range(len(c_op_list)):
		#print i
		c_op_list[i] = c_op_list[i].full()
		c_op_list[i]= scipy.sparse.csc_matrix(c_op_list[i])



	return c_op_list

#Defines the expectation operators that you want measured
def expectationOperators(lattice, settings):
	#Constructs the operators
	si = qutip.qeye(2)
	#sx = qutip.sigmax()
	#sy = qutip.sigmay()
	#sz = qutip.sigmaz()
	sp = qutip.sigmap()
	sm = qutip.sigmam()
	no = sp*sm
	#print no

	#Constructs operators for the chain of length N
	#si_list = []
	#sx_list = []
	#sy_list = []
	#sz_list = []
	#sp_list = []
	#sm_list = []
	no_list = []

	#si_list.append(qutip.tensor([si] * N))

	for k in range(lattice.numberOfSites):
		op_list = [si] * lattice.numberOfSites
		#op_list[k] = sx
		#sx_list.append(qutip.tensor(op_list))
		#op_list[k] = sy
		#sy_list.append(qutip.tensor(op_list))
		#op_list[k] = sz
		#sz_list.append(qutip.tensor(op_list))
		#op_list[k] = sp
		#sp_list.append(qutip.tensor(op_list))
		#op_list[k] = sm
		#sm_list.append(qutip.tensor(op_list))
		op_list[k] = no
		no_list.append(qutip.tensor(op_list))

	#Defines the expectation operators
	e_op_list = []

	#This adds on the measurement of the average population per site
	e_op_list.append(no_list[0])
	for i in range(1,lattice.numberOfSites):
		e_op_list[0] += no_list[i]
	e_op_list[0] /= lattice.numberOfSites

	for i in range(len(e_op_list)):
		e_op_list[i] = e_op_list[i].full()
		e_op_list[i] = scipy.sparse.csc_matrix(e_op_list[i])

	return e_op_list

def main(argv):
	#Sets up the lattice (add any command line inputs as shown)
	lattice = Lattice()
	lattice.numberOfSites = 5 #int(argv[0])
	#Sets the parameters of the Hamiltonian and the rates of the jump operators
	#This is not included in the backend to allow for as many rates as you like
	parameters = Parameters()
	parameters.omega = 1.0# float(argv[1])
	#Simulation settings
	settings = QJMCAA.Settings()
	settings.numberOfTrajectories = 1000
	approxAccuracy = 6
	dt = (settings.T/(settings.numberOfPoints-1))
	settings.smallestDt = dt*pow(10.0,-approxAccuracy)

	#Sets up how the data will be saved
	savingSettings = QJMCAA.SavingSettings()
	savingSettings.model = 'example'
	savingSettings.addParameter('Omega',parameters.omega)
	savingSettings.addExpectionSave('Population')

	#Gets the defined Hamiltonian
	H = Hamiltonian(parameters, lattice)
	#Gets the defined jump operators
	jumpOps = jumpOperators(parameters, lattice)
	#Gets the defined expectation operators
	eOps = expectationOperators(lattice, settings)

	#Gets the initial state you want to run it for
	psi0 = initialStateDefine(lattice.numberOfSites)

	#Runs the simulation
	QJMCAA.QJMCRun(settings, savingSettings, H, jumpOps, eOps, psi0)

if __name__ == "__main__":
	main(sys.argv[1:])
