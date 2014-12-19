#from gaussian_fitting import *
import numpy.random as R

mu0 = 0
sigma0 = 1
print("Real values: %s " % (mu0, sigma0))

## Compare the two methods
n = 1000
Mu1, Sig1 = zeros(n), zeros(n)
Mu2, Sig2 = zeros(n), zeros(n)
Mu3, Sig3 = zeros(n), zeros(n)
for i in range(n):
    s = R.normal(size=2000, loc=mu0, scale=sigma0)
    s = s[s>-0.8]
    #print "Empirical mean and std", s.mean(), s.std()
    m1,s1 = GF.gaussian_fit_cdf(s)
    m2,s2 = GF.gaussian_fit_hist(s, bins=10)
    m3,s3 = s.mean(), s.std()

    Mu1[i] = m1; Sig1[i] = s1
    Mu2[i] = m2; Sig2[i] = s2
    Mu3[i] = m3; Sig3[i] = s3

figure(1); title("Mu")
b = r_[mu0-1:mu0+1:40j]
kw = dict(alpha=0.4, histtype='stepfilled')
hist(Mu1, bins=b, label='PDF', **kw)
hist(Mu2, bins=b, label='Histo', **kw)
#hist(Mu3, bins=b, label='Moments', **kw)
legend(loc='best')
axvline(mu0, color='k', lw=2)
grid(True)

figure(2); title("Sigma")
b = r_[sigma0-1:sigma0+1:40j]
kw = dict(alpha=0.4, histtype='stepfilled')
hist(Sig1, bins=b, label='PDF', **kw)
hist(Sig2, bins=b, label='Histo', **kw)
#hist(Sig3, bins=b, label='Moments', **kw)
legend(loc='best')
axvline(sigma0, color='k', lw=2)
grid(True)

show()


