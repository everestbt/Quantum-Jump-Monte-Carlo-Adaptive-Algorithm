import QJMCSetUp

def evolvePsi(psi, t, dt, HEff):
    HEffExponentDt = QJMCSetUp.HEffExponentProduction(HEff,dt)
    psi = HEffExponentDt.dot(psi)
    t += dt
    return t, psi

def evolvePsiByExponent(psi, t, dt, HEffExponentDt):
    psi = HEffExponentDt.dot(psi)
    t += dt
    return t, psi
