#Edited 12/3/17 Ben Everest
#The file for the Quantum Jump Monte Carlo adapative algorithm.
#It contains the base classes which are used by the user to apply settings
#and it contains the code to run the algorithm
import numpy
import os
import progressbar
import random
import time
from datetime import timedelta

import QJMCMeasure
import QJMCSetUp
import QJMCMath
import QJMCJump
import QJMCEvolve
import QJMCJumpLowMemory

class Settings:
	def __init__(self):
		self.numberOfPoints = 1001
		self.T = 1.0
		self.numberOfTrajectories = 1
		self.location = os.path.dirname(os.path.realpath(__file__)) + '/'
		self.randomInitialStates = False
		self.smallestDt = 0.01
		#Not currently used
		self.memorySet = False
		self.memoryAvailable = 1
	def setMemory(val,measure):
		self.memoryAvailable=True
		if measure == 'B':
			self.memoryAvailable = val
		elif measure == 'KB':
			self.memoryAvailable = val * pow(10,3)
		elif measure == 'MB':
			self.memoryAvailable = val * pow(10,6)
		elif measure == 'GB':
			self.memoryAvailable = val * pow(10,9)

class SavingHistogram:
	def __init__(self):
		self.minimum = 0.0
		self.maximum = 1.0
		self.numberOfBars = 100

class SavingParameter:
	def __init__(self):
		self.name = ""
		self.value = 0.0

class SavingExpectation:
	def __init__(self):
		self.name = ''

class SavingSettings:
	def __init__(self):
		self.modelName = ''
		self.savingParameters = []
		self.expectationSave = []
		self.dataDecimalPlaces = 6
		self.histograms = False
		self.histogramOptions = []
	def addParameter(self,name,value):
		item = SavingParameter()
		item.name = name
		item.value = value
		self.savingParameters.append(item)
	def addExpectionSave(self,name):
		expect = SavingExpectation()
		expect.name = name
		self.expectationSave.append(expect)
	def addHistogram(self,minimum,maximum,numberOfBars):
		hist = SavingHistogram()
		hist.minimum = minimum
		hist.maximum = maximum
		hist.numberOfBars = numberOfBars
		self.histogramOptions.append(hist)

def QJMCRun(settings, savingSettings, H, jumpOps, eOps, psi0):
	#Tests the inputs
	QJMCSetUp.dimensionTest(H,jumpOps,eOps,psi0)
	QJMCSetUp.typeTest(settings, savingSettings, H, jumpOps, eOps, psi0)

	startTime = time.time()

	#Creates the time array
	tList = numpy.linspace(0,settings.T,settings.numberOfPoints)

	#Gets the pairs of jump operators for later use
	jumpOpsPaired = QJMCSetUp.jumpOperatorsPaired(jumpOps)

	#Produces expectation operators squared to help produce variance
	QJMCSetUp.addExpectationSquared(eOps)

	#Produces results arrays
	eResults = QJMCSetUp.eResultsProduction(eOps,settings.numberOfPoints)

	#Produces the histograms list
	if (savingSettings.histograms):
		histograms = QJMCSetUp.histogramProduction(savingSettings.histogramOptions,
			settings.numberOfPoints)
	else:
		histograms = []

	#Produces the effective hamiltonain as an exponent and gets the smaller set
	HEffExponentDtSet, dtSet = QJMCSetUp.HEffExponentSetProductionBinary(H,
		jumpOpsPaired, tList[1], settings)

	#Used to track the progress of the simulation
	bar = progressbar.ProgressBar(maxval=settings.numberOfTrajectories, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for traj in range(settings.numberOfTrajectories):
		bar.update(traj)
		#Resets the system
		#Selects a random initial state if requested
		if settings.randomInitialStates:
			psi0 = QJMCSetUp.randomInitialState(H)
		psi = psi0
		previousPsi = psi
		t = 0
		index = 0

		#Performs the first measure at t=0
		index = QJMCMeasure.measure(index, psi, eResults, eOps,
			histograms, savingSettings)
		#Selects the random number for the first time selection
		r = random.random()
		#START OF LOOP
		while (index < settings.numberOfPoints):
			#Saves the previous state if necessary
			previousPsi = psi
			#Progresses to the next time step
			#TODO change to evolve
			psi = HEffExponentDtSet[0].dot(psi)
			#Calculates the survivial probability
			survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)
			#Checks if it has jumped in the last time step
			if survivialProbability < r:
				#Performs the jump at the more accurate time
				t, psi, r = QJMCJump.jumpBinary(tList[index - 1], previousPsi,
					r, jumpOps, jumpOpsPaired, dtSet, HEffExponentDtSet)
			#Performs the measurement
			index = QJMCMeasure.measure(index, psi, eResults,eOps,
				histograms, savingSettings)
	bar.update(settings.numberOfTrajectories)
	print('Saving')
	#Averages all the results
	eResults = QJMCMeasure.averageResults(eResults,settings)
	if (savingSettings.histograms):
		histograms = QJMCMeasure.averageHistograms(histograms,settings.numberOfTrajectories)
	#Saves the results
	QJMCMeasure.saveResults(settings,savingSettings,eResults)
	if (savingSettings.histograms):
		QJMCMeasure.saveHistograms(settings,savingSettings,histograms)

	endTime = time.time()
	print('Time taken ' + str(timedelta(seconds=endTime-startTime)))

def QJMCRunLowMemory(settings, savingSettings, tList, H, jumpOps, eOps, psi0):
	#Tests the inputs
	QJMCSetUp.dimensionTest(H,jumpOps,eOps,psi0)
	#TODO add type test on tList
	QJMCSetUp.typeTest(settings, savingSettings, H, jumpOps, eOps, psi0)

	startTime = time.time()

	#Sets the number of points based on the tList
	#TODO check this works correctly
	settings.numberOfPoints = len(tList)

	#Gets the pairs of jump operators for later use
	jumpOpsPaired = QJMCSetUp.jumpOperatorsPaired(jumpOps)

	#Gets the effective Hamiltonian
	HEff = QJMCSetUp.HEffProduction(H, jumpOpsPaired)

	#Produces expectation operators squared to help produce variance
	QJMCSetUp.addExpectationSquared(eOps)

	#Produces results arrays
	eResults = QJMCSetUp.eResultsProduction(eOps,settings.numberOfPoints)

	#Produces the histograms
	if (savingSettings.histograms):
		histograms = QJMCSetUp.histogramProduction(savingSettings.histogramOptions,
			settings.numberOfPoints)

	#Used to track the progress of the simulation
	bar = progressbar.ProgressBar(maxval=settings.numberOfTrajectories, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for traj in range(settings.numberOfTrajectories):
		bar.update(traj)
		#Resets the system
		#Selects a random initial state if requested
		if settings.randomInitialStates:
			psi0 = QJMCSetUp.randomInitialState(H)
		psi = psi0
		t = 0
		index = 0

		#Performs the first measure at t=0
		index = QJMCMeasure.measure(index, psi, eResults, eOps,
			histograms, savingSettings)
		#Selects the random number for the first time selection
		r = random.random()
		#START OF LOOP
		while (index < settings.numberOfPoints):
			#Saves the previous state if necessary
			previousPsi = psi
			#Progresses to the next time step
			dt = tList[index] - tList[index-1]
			t, psi = QJMCEvolve.evolvePsi(psi, t, dt, HEff)
			#Calculates the survivial probability
			survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)
			#Checks if it has jumped in the last time step
			if survivialProbability < r:
				#Performs the jump at the more accurate time
				t, psi = QJMCJumpLowMemory.jump(tList[index - 1], dt,
					settings.smallestDt, psi,r, jumpOps, jumpOpsPaired, HEff)
				#Catches up the accurate time to the index time, performing any
				#other jumps
				t, psi, r = QJMCJumpLowMemory.catchUpApprox(t, tList[index],
					settings.smallestDt, psi, jumpOps, jumpOpsPaired,HEff)
			#Performs the measurement
			index = QJMCMeasure.measure(index, psi, eResults,eOps,
				histograms, savingSettings)
	bar.update(settings.numberOfTrajectories)
	print('Saving')
	#Averages all the results
	eResults = QJMCMeasure.averageResults(eResults,settings)
	if (savingSettings.histograms):
		histograms = QJMCMeasure.averageHistograms(histograms,settings.numberOfTrajectories)
	#Saves the results
	QJMCMeasure.saveResults(settings,savingSettings,eResults)
	if (savingSettings.histograms):
		QJMCMeasure.saveHistograms(settings,savingSettings,histograms)

	endTime = time.time()
	print('Time taken ' + str(timedelta(seconds=endTime-startTime)))
