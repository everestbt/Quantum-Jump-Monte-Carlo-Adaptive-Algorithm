import QJMCSetUp

#TODO use this where possible in the low memory solution
def evolvePsi(psi, t, dt, HEff):
    HEffExponentDt = QJMCSetUp.HEffExponentProduction(HEff,dt)
    psi = HEffExponentDt.dot(psi)
    t += dt
    return t, psi
