#Edited 12/3/17 Ben Everest
#Unit tests for the set-up functions
import qutip
import numpy
import scipy
import unittest

import QJMCSetUp
import QJMCAA
import QJMCMath


class TestQJMCSetUp(unittest.TestCase):
    def test_jumpOperatorsPaired(self):
        #Sets up the jump operators
        sp = qutip.sigmap()
    	sm = qutip.sigmam()
    	no = sp*sm

    	jumpOps = []
        jumpOps.append(sm)
        jumpOps.append(sp)
        jumpOps.append(no)

        for i in range(len(jumpOps)):
    		jumpOps[i] = jumpOps[i].full()
    		jumpOps[i]= scipy.sparse.csc_matrix(jumpOps[i])
        #Runs the function
        jumpOpsPaired = QJMCSetUp.jumpOperatorsPaired(jumpOps)

        #Calculating what they should be
        jumpOpsExpect = []

        jumpOpsExpect.append(scipy.sparse.csc_matrix([[complex(1,0),0.0],[0.0,0.0]]))
        jumpOpsExpect.append(scipy.sparse.csc_matrix([[0.0,0.0],[0.0,complex(1,0)]]))
        jumpOpsExpect.append(scipy.sparse.csc_matrix([[complex(1,0),0.0],[0.0,0.0]]))

        fail = -1
        for i in range(len(jumpOps)):
            if (jumpOpsPaired[i] - jumpOpsExpect[i]).nnz == complex(0,0):
                continue
            else:
                fail = i
                break
        self.assertEqual(fail,-1)

    def test_typeTest(self):
        #Test the settings
        settings = QJMCAA.Settings()
        savingSettings = QJMCAA.SavingSettings()

        H = scipy.sparse.csc_matrix((3, 4), dtype=numpy.int8)

        jumpOps = []
        jumpOps.append(scipy.sparse.csc_matrix((3, 4), dtype=numpy.int8))

        eOps = []
        eOps.append(scipy.sparse.csc_matrix((3, 4), dtype=numpy.int8))

        psi0 = numpy.ndarray(shape=(2,2), dtype=float, order='F')

        QJMCSetUp.typeTest(settings,savingSettings,H,jumpOps,eOps,psi0)

        self.assertTrue(True)

    def test_randomInitialState(self):
        sx = qutip.sigmax()

    	H =  sx

    	#All objects must be in a scipy sparse format
    	H = H.full()
    	H = scipy.sparse.csc_matrix(H)

        psi0 = []
        num = 100
        for i in range(num):
            psi0.append(QJMCSetUp.randomInitialState(H))

        for i in range(num):
            self.assertAlmostEqual(1.0,QJMCMath.calculateSquareOfWavefunction(psi0[i]))

        for i in range(num - 1):
            for j in range(H.get_shape()[0]):
                self.assertNotAlmostEqual(psi0[i][j],psi0[i+1][j])

    def test_dimensionTest(self):
        H = scipy.sparse.csc_matrix((2, 2), dtype=numpy.int8)

        jumpOps = []
        jumpOps.append(scipy.sparse.csc_matrix((2, 2), dtype=numpy.int8))

        eOps = []
        eOps.append(scipy.sparse.csc_matrix((2, 2), dtype=numpy.int8))

        psi0 = numpy.ndarray(shape=(2,1), dtype=float, order='F')

        QJMCSetUp.dimensionTest(H,jumpOps,eOps,psi0)
        self.assertTrue(True)

if __name__ == '__main__':
	unittest.main()
