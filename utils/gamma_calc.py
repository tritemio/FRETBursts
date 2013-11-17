#
# 1/S = OME + SIG*Epr  (Epr = Fad/(Fad+Fdd))
#
#
SL = 0.5
SH = 0.5
EL = 0.21
EH = 0.75

DeltaS = 0.02

def cal_gamma(SL,SH,EL,EH):
    SIG = (1./SH - 1./SL)/(EH-EL)
    OME = 1./SL - SIG*EL
    gamma = (OME-1)/(OME+SIG-1)
    return gamma

def cal_gamma_beta(SL,SH,EL,EH):
    SIG = (1/SH - 1/SL)/(EH-EL)
    OME = 1/SL - SIG*EL
    gamma = (OME-1)/(OME+SIG-1)
    beta = OME+SIG-1
    return gamma, beta


SL1, SL2 = SL-DeltaS, SL+DeltaS
SH1, SH2 = SH-DeltaS, SH+DeltaS

g = cal_gamma(SL,SH,EL,EH)
g1 = cal_gamma(SL2,SH1,EL,EH)
g2 = cal_gamma(SL1,SH2,EL,EH)

print g,g1,g2
#print b,b1,b2

f = lambda DeltaS, S, EL, EH: cal_gamma(S+DeltaS, S-DeltaS, EL, EH)

plot(d, f(d,0.7,0.2,0.2+0.5))
