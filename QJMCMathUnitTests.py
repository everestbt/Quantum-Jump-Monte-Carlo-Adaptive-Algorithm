#Edited 12/3/17 Ben Everest
#Unit tests for the math functions
import numpy
import unittest
import math

import QJMCMath

class TestQJMCMath(unittest.TestCase):

    def test_calculateSquareOfWavefunction(self):
        psi = numpy.ndarray(shape=(2,1),dtype = complex)

        psi[0,0]=complex(1,0)
        psi[1,0] = complex(0,0)
        P = QJMCMath.calculateSquareOfWavefunction(psi)
        self.assertAlmostEqual(P,1.0)

        psi[0,0]=complex(0.5,0)
        psi[1,0] = complex(0,0)
        P = QJMCMath.calculateSquareOfWavefunction(psi)
        self.assertAlmostEqual(P,0.25)

        psi[0,0]=complex(0.5,0)
        psi[1,0] = complex(0.5,0)
        P = QJMCMath.calculateSquareOfWavefunction(psi)
        self.assertAlmostEqual(P,0.5)

        psi[0,0]=complex(0.5,0)
        psi[1,0] = complex(0,0.5)
        P = QJMCMath.calculateSquareOfWavefunction(psi)
        self.assertAlmostEqual(P,0.5)


    def test_neighbour(self):
        N = 10
        #Zero shift
        pos = QJMCMath.neighbour(1, N, 0)
        self.assertEqual(pos,1)
        #1 shift
        pos = QJMCMath.neighbour(1, N, 1)
        self.assertEqual (pos,2)
        #-1 shift
        pos = QJMCMath.neighbour(1, N, -1)
        self.assertEqual (pos,0)
        #2 shift
        pos = QJMCMath.neighbour(1, N, 2)
        self.assertEqual (pos,3)
        #-2 shift
        pos = QJMCMath.neighbour(1, N, -2)
        self.assertEqual (pos,(N-1))
        #N shift
        pos = QJMCMath.neighbour(1, N, N)
        self.assertEqual (pos,1)
        #-N shift
        pos = QJMCMath.neighbour(1, N, -N)
        self.assertEqual (pos,1)

    def test_normalise(self):
        psi = [[complex(1,0)],[complex(1,0)]]
        psiN = QJMCMath.normalise(psi)
        norm = math.sqrt(QJMCMath.calculateSquareOfWavefunction(psiN))
        self.assertAlmostEqual(norm,1.0)
        psiTest = [[complex(1.0/math.sqrt(2.0),0)],[complex(1.0/math.sqrt(2.0),0)]]
        self.assertItemsEqual(psiN,psiTest)

        psi = [[complex(1,1)],[complex(1,1)]]
        psiN = QJMCMath.normalise(psi)
        norm = math.sqrt(QJMCMath.calculateSquareOfWavefunction(psiN))
        self.assertAlmostEqual(norm,1.0)
        psiTest = [[complex(0.5,0.5)],[complex(0.5,0.5)]]
        self.assertItemsEqual(psiN,psiTest)

        psi = [[complex(1,0)],[complex(0,0)]]
        psiN = QJMCMath.normalise(psi)
        norm = math.sqrt(QJMCMath.calculateSquareOfWavefunction(psiN))
        self.assertAlmostEqual(norm,1.0)
        psiTest = [[complex(1.0,0)],[complex(0,0)]]
        self.assertItemsEqual(psiN,psiTest)

if __name__ == '__main__':
    unittest.main()
