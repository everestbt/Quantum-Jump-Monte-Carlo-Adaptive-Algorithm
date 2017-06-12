import qutip
import unittest
import numpy
import scipy

import QJMCMeasure
import QJMCAA

class TestQJMCMeasure(unittest.TestCase):
    def test_measure(self):
        #Defines all needed variables for a test
        index = 0

        psi = qutip.basis(2, 1)
    	psi = psi.full()

        eResults = []

        eResults.append(numpy.zeros(1))

        sp = qutip.sigmap()
    	sm = qutip.sigmam()
    	no = sp*sm
        eOps = []
    	eOps.append(no)
        for i in range(len(eOps)):
    		eOps[i] = eOps[i].full()
    		eOps[i] = scipy.sparse.csc_matrix(eOps[i])

        histograms = []

        savingSettings = QJMCAA.SavingSettings()

        #Tests with no histograms set
        QJMCMeasure.measure(index, psi, eResults, eOps, histograms, savingSettings)

        self.assertEqual(histograms,[])
        self.assertEqual(eResults,[[0.]])
