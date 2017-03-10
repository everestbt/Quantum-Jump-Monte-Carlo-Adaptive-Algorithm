import numpy
import QJMCMath

#Performs all the measurements
def measure(index, psi, eResults,eOps, histograms, savingSettings):
	#Extracts normalised wavefunctions
	psiNormalised = QJMCMath.normalise(psi)
	psiNormalisedLeft = numpy.conjugate(numpy.transpose(psiNormalised))[0]

	#Does the measurements:
	for i in range(len(eResults)):
		value = numpy.asscalar(numpy.real(numpy.dot(psiNormalisedLeft,
			(eOps[i].dot(psiNormalised)))))
		eResults[i][index] += value
		#Checks whether histograms are wanted
		if (savingSettings.histograms):
			#Goes through all histograms requested
			if (i < len(savingSettings.histogramOptions)):
				#Extracts the values for that histogram
				minimum = savingSettings.histogramOptions[i].minimum
				maximum = savingSettings.histogramOptions[i].maximum
				numberOfBars = savingSettings.histogramOptions[i].numberOfBars
				#Gets the position in the histogram
				position = int(((value - minimum)/maximum) * numberOfBars)
				#If the value lies outside of the range, places it on the edges
				if (position >= numberOfBars):
					position = numberOfBars - 1
				elif (position < 0):
					position = 0
				#Puts the count into the correct histogram
				histograms[i][index][position] += 1

	#Advances the measurement by 1
	index += 1
	return index

def averageResults(eResults,settings):
	for i in range(len(eResults)):
		eResults[i] = eResults[i]/settings.numberOfTrajectories

	return eResults

def averageHistograms(histograms,numberOfTrajectories):
	for i in range(len(histograms)):
		histograms[i]  = histograms[i] / numberOfTrajectories
	return histograms

def variance(av, avSquared):
	return (avSquared - numpy.power(av,2))

def nameTheFile(savingSettings, i):
	name = savingSettings.modelName
	for item in savingSettings.savingParameters:
		name += item.name
		name +=str(item.value)
	name += savingSettings.expectationSave[i].name
	name += '.txt'
	return name

def saveResults(settings,savingSettings,eResults):
	#The start of the squared results
	sqStart = len(eResults)/2
	#Produces the time list
	tList = numpy.linspace(0,settings.T,settings.numberOfPoints)
	#Saves each one that was requested by the user
	for i in range(len(savingSettings.expectationSave)):
		name = nameTheFile(savingSettings,i)
		var = variance(eResults[i],eResults[i+sqStart])
		numpy.savetxt(name, numpy.c_[tList,eResults[i],var],
			fmt='%.'+str(savingSettings.dataDecimalPlaces)+'e')

def saveHistograms(settings,savingSettings,histograms):
	#Produces the time list
	tList = numpy.linspace(0,settings.T,settings.numberOfPoints)
	for i in range(len(histograms)):
		name = nameTheFile(savingSettings,i)
		#Corrects the name to include Histogram
		name = name[:-4]
		name += 'Histogram'
		name += '.txt'
		open(name, 'w').close()
		with open(name,'a') as f:
			numpy.savetxt(f, histograms[i],
				fmt='%.'+str(savingSettings.dataDecimalPlaces)+'e')
