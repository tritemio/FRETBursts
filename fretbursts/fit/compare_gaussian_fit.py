from gaussian_fitting import *

mu0 = -0.1
sigma0 = 1.8
print("Real values: %s %s" % (mu0, sigma0))

## Compare the two methods
n = 10000
Mu, Sigma = zeros(n), zeros(n)
Muh, Sigmah = zeros(n), zeros(n)
Mu1, Sigma1 = zeros(n), zeros(n)
for i in range(n):
    s = R.normal(size=100, loc=mu0, scale=sigma0)
    #print "Empirical mean and std", s.mean(), s.std()
    mu, sigma =  gaussian_fit(s, [0,1])
    muh, sigmah, h, bins =  gaussian_fit_hist(s, [0,1])
    mu1, sigma1 = s.mean(), s.std()

    Mu[i] = mu; Sigma[i] = sigma
    Muh[i] = muh; Sigmah[i] = sigmah
    Mu1[i] = mu1; Sigma1[i] = sigma1

figure(1); title("Mu")
b = r_[mu0-1:mu0+1:40j]
hist(Mu, bins=b, alpha=0.3, label='ECDF')
hist(Muh, bins=b, alpha=0.3, color='r', label='Histo')
hist(Mu1, bins=b, alpha=0.3, color='y', label='Moments')
legend(loc='best')
axvline(mu0, color='k', lw=2)
grid(True)

figure(2); title("Sigma")
b = r_[sigma0-1:sigma0+1:40j]
hist(Sigma, bins=b, alpha=0.3, label='ECDF')
hist(Sigmah, bins=b, alpha=0.3, color='r', label='Histo')
hist(Sigma1, bins=b, alpha=0.3, color='y', label='Moments')
legend(loc='best')
axvline(sigma0, color='k', lw=2)
grid(True)

show()


