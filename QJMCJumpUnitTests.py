#Edited 12/3/17 Ben Everest
#Unit tests for the set-up functions
import qutip
import scipy
import unittest

import QJMCSetUp
import QJMCAA
import QJMCJump
import QJMCEvolve


class TestQJMCJump(unittest.TestCase):
	def test_jumpBinary(self):
		tStart = 0
		psi0 = qutip.basis(2, 0)
		psi0 = psi0.full()

		sp = qutip.sigmap()
		sm = qutip.sigmam()
		jumpOps = []
		#Decay
		jumpOps.append(5 * sm)
		jumpOps.append(5 * sp)
		for i in range(len(jumpOps)):
			jumpOps[i] = jumpOps[i].full()
			jumpOps[i]= scipy.sparse.csc_matrix(jumpOps[i])
		jumpOpsPaired = QJMCSetUp.jumpOperatorsPaired(jumpOps)

		sx = qutip.sigmax()
		H = sx
		#All objects must be in a scipy sparse format
		H = H.full()
		H = scipy.sparse.csc_matrix(H)

		settings = QJMCAA.Settings()
		settings.smallestDt = 0.0001

		HEffExponentDtSet, dtSet = QJMCSetUp.HEffExponentSetProductionBinary(H,
			jumpOpsPaired, 0.1, settings)

		#print(dtSet)
		t, psi = QJMCEvolve.evolvePsiByExponent(psi0, tStart, dtSet[0], HEffExponentDtSet[0])

		for _ in range(1000):
			r = 0.5
			t, _, r = QJMCJump.jumpBinary(tStart, psi0, r, jumpOps, jumpOpsPaired,
		 		dtSet, HEffExponentDtSet)
			self.assertAlmostEqual(t,0.1)


if __name__ == '__main__':
	unittest.main()
