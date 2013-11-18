"""
Simulate FRET distribution using synthetic bursts.

NOTE: only ONE burst size is considered.
"""
from scipy.stats import poisson
from scipy.misc import factorial

lam = 3
SSS = 50
Sd = SSS
Sa = SSS
Pa = poisson(lam)

NSamples = 1e5
BG = Pa.rvs(NSamples)
dBG = BG-lam

## Real extension of poisson PMF
poiss_pdf = lambda x,l: (exp(-l)*(l**(x)))/(factorial(x))
## Real extension of poisson PMF with mean subtraction
poiss_pdf_m = lambda x,l: (exp(-l)*(l**(x+l)))/(factorial(x+l))

# Sample a reasonable (dense) range for poisson dist.
# step is 1/(2**n) to sample the integers
x = r_[:lam+5*sqrt(lam)+1:0.125*0.5] 
x_m = x-lam # axis for mean-substracted Poisson

plot(x_m,Pa.pmf(x_m+lam), lw=2)
hist(dBG, r_[-10:20]-0.5, histtype='step', normed=True, lw=2)
plot(x_m,poiss_pdf_m(x_m,lam), lw=2)
legend(["PMF (integers)", "Real extension of PMF", "Hist of 1e5 samples"], 
        loc='best')
grid(True)
kjhafdk
## D-only FRET 
E = 1.*dBG/(dBG + Sd) # MC
Ef = 1.*x_m/(x_m+Sd)      # Analytical

step = 0.0025
bins = r_[-1:2:step]+(0.5*step)
figure(figsize=(8,4))
hist(E, bins, histtype='bar', normed=True, alpha=0.3)
plot(Ef, poiss_pdf_m(x_m,lam)/step, lw=1.5)
xlim(-0.2,0.2)
if Sd/lam < 7: xlim(-0.4,0.4)
xlabel('E'); ylabel("Burst PDF"); grid(True)
title("E=0: (Burst PDF)*step [step=%.4f], Sd=%d, l_bga=%.1f" % (step, Sd, lam))
legend(["Analitycal","Simulation with %.e bursts" %NSamples], loc='best')
#savefig("E=0: FRET sim MC vs Analitycal, Sd=%d, BGa=%d.png" %(Sd, lam))
kjfhsd
## 100% FRET (here dBG represent BG in D-ch)
E1 = 1.*Sa/(Sa + dBG) # MC
Ef1 = 1.*Sa/(Sa+x_m)    # Analytical

step = 0.0025
bins1 = r_[-1:2:step]+(0.5*step)
figure(figsize=(8,4))
hist(E1, bins1, histtype='bar', normed=True, alpha=0.3)
plot(Ef1, poiss_pdf_m(x_m,lam)/step, lw=1.5)
xlim(0.8,1.2)
if Sa/lam < 7: xlim(0.6,1.4)
xlabel('E'); ylabel("Burst PDF"); grid(True)
title("E=1: (Burst PDF)*step [step=%.4f], Sa=%d, l_bgd=%.1f" % (step, Sd, lam))
legend(["Analitycal","Simulation with %.e bursts" %NSamples], loc='best')
savefig("E=1: FRET sim MC vs Analitycal, Sa=%d, BGd=%d.png" %(Sa, lam))


#from glob import glob
#from subprocess import call
#l = glob("E=*=5.png")
#l = [l[i] for i in [0,3,2,5,1,4]]
#call(["montage"]+l+"-tile 2x3 -geometry 800x400+1+1 mont_BG=5.png".split(" "))
