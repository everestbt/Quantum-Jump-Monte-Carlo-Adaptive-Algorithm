import numpy
import scipy
import os
import progressbar
import random
import time
import sys
from datetime import timedelta

import QJMCMeasure
import QJMCSetUp
import QJMCMath
import QJMCJump

class Settings:
	def __init__(self):
		self.numberOfPoints = 1001
		self.T = 1.0
		self.numberOfTrajectories = 1
		self.location = os.path.dirname(os.path.realpath(__file__)) + '/'
		self.accuracyMagnitude = 3
		self.randomInitialStates = False

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

	#Produces the histograms
	if (savingSettings.histograms):
		histograms = QJMCSetUp.histogramProduction(savingSettings.histogramOptions,
			settings.numberOfPoints)

	#Produces the effective hamiltonain as an exponent and gets the smaller set
	HEffExponentDt, HEffExponentDtSet = QJMCSetUp.HEffExponentProduction(H,
		jumpOpsPaired, tList[1], settings.accuracyMagnitude)

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
			psi = HEffExponentDt.dot(psi)
			#Progresses the time, is commented out as it is not needed since
			#the index tracks the actual time
			#t+=tList[1]
			#Calculates the survivial probability
			survivialProbability = QJMCMath.calculateSquareOfWavefunction(psi)
			#Checks if it has jumped in the last time step
			if survivialProbability < r:
				#Performs the jump at the more accurate time
				t, psi = QJMCJump.jump(tList[index - 1], tList[1], previousPsi,
					r, jumpOps, jumpOpsPaired, HEffExponentDtSet)
				#Catches up the accurate time to the index time, performing any
				#other jumps
				t, psi, r = QJMCJump.catchUpApprox(t, tList[index], tList[1], psi,
					jumpOps, jumpOpsPaired,HEffExponentDtSet)
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