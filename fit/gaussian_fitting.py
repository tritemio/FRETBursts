"""
This module provides functions to fit gaussian distributions and gaussian 
distribution mixtures (2 components). Several fittings methods are provided.

Single Gaussian distribution fit:
    gaussian_fit_hist()
    gaussian_fit_cdf()
    gaussian_fit_pdf()

Mixture of 2 Gaussian distribution fit:
    two_gaussian_fit_hist()
    two_gaussian_fit_cdf()
    two_gaussian_fit_EM()

Also some functions to fit 2-D gaussian distributions and mixtures are
implemented (needs more testing).
"""

from pylab import normpdf
from numpy import *
import numpy.random as R
import scipy.optimize as O
import scipy.stats as S

from scipy.special import erf
from scipy.optimize import leastsq, minimize
from scipy.ndimage.filters import gaussian_filter1d

##
# Single gaussian distribution
#

def gaussian_fit_curve(x, y, mu0=0, sigma0=1, a0=None, return_all=False, 
        **kwargs):
    """Gaussian fit of curve (x,y).
    If a0 is None then only (mu,sigma) are fitted (to a gaussian density).
    `kwargs` are passed to the leastsq() function.

    If return_all=False then return only the fitted (mu,sigma) values
    If return_all=True (or full_output=True is passed to leastsq)
    then the full output of leastsq is returned.
    """
    if a0 is None:
        gauss_pdf = lambda x,m,s: exp(-((x-m)**2)/(2*s**2))/(sqrt(2*pi)*s)
        err_fun = lambda p,x,y: gauss_pdf(x,*p) - y
        res = leastsq(err_fun, x0=[mu0,sigma0], args=(x,y), **kwargs)
    else:
        gauss_fun = lambda x,m,s,a: a*sign(s)*exp(-((x-m)**2)/(2*s**2))
        err_fun = lambda p,x,y: gauss_fun(x,*p) - y
        res = leastsq(err_fun, x0=[mu0,sigma0,a0], args=(x,y), **kwargs)

    if 'full_output' in kwargs: return_all = kwargs['full_output']
    mu, sigma = res[0][0], res[0][1]
    if return_all: return res
    return mu, sigma 

def get_epdf(s, smooth=0, N=1000, smooth_pdf=False, smooth_cdf=True):
    """Compute the empirical PDF of the sample `s`.

    If smooth > 0 then apply a gaussian filter with sigma=smooth.
    N is the number of points for interpolation of the CDF on a uniform range.
    """
    ecdf = [sort(s), arange(0.5,s.size+0.5)*1./s.size]
    #ecdf = [sort(s), arange(s.size)*1./s.size]
    _x = linspace(s.min(),s.max(),N)
    ecdfi = [_x, interp(_x, ecdf[0],ecdf[1])]
    if smooth_cdf and smooth > 0:
        ecdfi[1] =gaussian_filter1d(ecdfi[1], sigma=smooth)
    epdf = [ecdfi[0][:-1], diff(ecdfi[1])/diff(ecdfi[0])]
    if smooth_pdf and smooth > 0:
        epdf[1] =gaussian_filter1d(epdf[1], sigma=smooth)
    return epdf

def gaussian_fit_pdf(s, mu0=0, sigma0=1, a0=1, return_all=False, 
        leastsq_kwargs={}, **kwargs):
    """Gaussian fit of samples s using a fit to the empirical PDF.
    If a0 is None then only (mu,sigma) are fitted (to a gaussian density).
    `kwargs` are passed to get_epdf().
    If return_all=False then return only the fitted (mu,sigma) values
    If return_all=True (or full_output=True is passed to leastsq)
    then the full output of leastsq and the PDF curve is returned.
    """
    ## Empirical PDF
    epdf = get_epdf(s, **kwargs) 
   
    res = gaussian_fit_curve(epdf[0], epdf[1], mu0, sigma0, a0, return_all,
            **leastsq_kwargs)
    if return_all: return res, epdf
    return res

def gaussian_fit_hist(s, mu0=0, sigma0=1, a0=None, bins=r_[-0.5:1.5:0.001],
        return_all=False, leastsq_kwargs={}, weights=None, **kwargs):
    """Gaussian fit of samples s fitting the hist to a Gaussian function.
    If a0 is None then only (mu,sigma) are fitted (to a gaussian density).
    kwargs are passed to the histogram function.
    If return_all=False then return only the fitted (mu,sigma) values
    If return_all=True (or full_output=True is passed to leastsq)
    then the full output of leastsq and the histogram is returned.
    `weights` optional weights for the histogram.
    """
    histogram_kwargs = dict(bins=bins, density=True, weights=weights)
    histogram_kwargs.update(**kwargs)
    H = histogram(s, **histogram_kwargs)
    x, y = 0.5*(H[1][:-1] + H[1][1:]), H[0]
    #bar(H[1][:-1], H[0], H[1][1]-H[1][0], alpha=0.3)
    
    res = gaussian_fit_curve(x, y, mu0, sigma0, a0, return_all,
            **leastsq_kwargs)
    if return_all: return res, H, x, y
    return res

def gaussian_fit_cdf(s, mu0=0, sigma0=1, return_all=False, **leastsq_kwargs):
    """Gaussian fit of samples s fitting the empirical CDF.
    Additional kwargs are passed to the leastsq() function.
    If return_all=False then return only the fitted (mu,sigma) values
    If return_all=True (or full_output=True is passed to leastsq)
    then the full output of leastsq and the histogram is returned.
    """
    ## Empirical CDF
    ecdf = [sort(s), arange(0.5,s.size+0.5)*1./s.size]
    
    ## Analytical Gaussian CDF
    gauss_cdf = lambda x, mu, sigma: 0.5*(1+erf((x-mu)/(sqrt(2)*sigma)))
    
    ## Fitting the empirical CDF
    err_func = lambda p, x, y: y - gauss_cdf(x, p[0], p[1])
    res = leastsq(err_func, x0=[mu0, sigma0], args=(ecdf[0],ecdf[1]),
            **leastsq_kwargs)
    if return_all: return res, ecdf
    return res[0]

def gaussian_fit_ml(s, mu_sigma_guess=[0.5,1]):
    """Gaussian fit of samples s using the Maximum Likelihood (ML method).
    Didactical, since S.norm.fit() implements the same method.    
    """
    n = s.size
    ## Log-likelihood (to be maximized)
    log_l = lambda mu, sig: -n/2.*log(sig**2) -1./(2*sig**2)*sum((s-mu)**2)
    
    ## Function to be minimized
    min_fun = lambda p: -log_l(p[0], p[1])     
    
    res = O.minimize(min_fun, [0,0.5], method='powell',
                     options={'xtol': 1e-6, 'disp': True, 'maxiter': 1e9})
                     
    print res 
    mu, sigma = res['x']
    return mu, sigma

##
# Two-component gaussian mixtures
#

def two_gauss_mix_pdf(x, p):
    """PDF for the distribution of a mixture of two Gaussians."""
    mu1, sig1, mu2, sig2, a = p
    return a*normpdf(x, mu1, sig1) + (1-a)*normpdf(x, mu2, sig2)

def reorder_parameters(p):
    """Reorder 2-gauss mix params to have the 1st component with smaller mean.
    """
    if p[0] > p[2]:
        p = p[array([2,3,0,1,4])] # swap (mu1, sig1) with (mu2, sig2)
        p[4] = 1 - p[4]           # "swap" the alpha of the mixture
    return p

def two_gaussian_fit_EM(s, p0=[0,0.1,0.6,0.1,0.5], max_iter=300, ptol=1e-4,
        fix_mu=[0,0], fix_sig=[0,0], debug=False, w=None):
    """
    Fit the sample s with two gaussians using Expectation Maximization.
    `ptol`: convergence condition. Relative max variation of any parameter.
    `max_iter`: max number of iteration in case of non convergence.
    `fix_mu`: allow to fix the mean (mu) of a component (1 or True to fix)
    `fix_sig`: allow to fix the std. dev of a component (1 or True to fix)
    `w`: optional weigths, same size as `s` (for ex. 1/sigma^2 ~ nt).
    """
    assert size(p0) == 5
    if w is None: w = ones(s.size)
    assert w.size == s.size
    w *= (1.*w.size)/w.sum() # Normalize to (#samples), not to 1
    #w /= w.sum() # Normalize to 1
    if debug: assert abs(w.sum() - s.size) < 1e-6
    
    # Initial guess of parameters and initializations
    mu = array([p0[0], p0[2]])
    sig = array([p0[1], p0[3]])
    pi_ = array([p0[4], 1-p0[4]])

    gamma = zeros((2, s.size))
    N_ = zeros(2)
    p_new = array(p0)

    # EM loop
    counter = 0
    stop_iter, converged = False, False
    while not stop_iter:
        # Compute the responsibility func. (gamma) and the new parameters
        for k in [0,1]:
            gamma[k,:] = w*pi_[k]*normpdf(s, mu[k], sig[k]) / \
                    two_gauss_mix_pdf(s, p_new)
            ## Uncomment for SCHEME2
            #gamma[k,:] = pi_[k]*normpdf(s, mu[k], sig[k]) / \
            #        two_gauss_mix_pdf(s, p_new)
            N_[k] = gamma[k,:].sum()
            if not fix_mu[k]: 
                mu[k] = sum(gamma[k]*s)/N_[k]
                ## Uncomment for SCHEME2
                #mu[k] = sum(w*gamma[k]*s)/N_[k]
            if not fix_sig[k]:
                sig[k] = sqrt( sum(gamma[k]*(s-mu[k])**2)/N_[k] )
            pi_[k] = N_[k]/s.size
        p_old = p_new
        p_new = array([mu[0], sig[0], mu[1], sig[1], pi_[0]])
        if debug:
            assert abs(N_.sum() - s.size)/float(s.size) < 1e-6 
            assert abs(pi_.sum() - 1) < 1e-6
        
        # Convergence check
        counter += 1
        fixed = concatenate([fix_mu, fix_sig, [0]]).astype(bool)
        relative_delta = abs(p_new[-fixed]-p_old[-fixed])/p_new[-fixed]
        converged = relative_delta.max() < ptol
        stop_iter = converged or (counter >= max_iter)
    
    if debug:
        print "Iterations: ", counter
    if not converged:
        print "WARNING: Not converged, max iteration (%d) reached." % max_iter
    return reorder_parameters(p_new)

def two_gaussian_fit_hist(s, bins=r_[-0.5:1.5:0.001], weights=None, 
        p0=[0.2,1,0.8,1,0.3], fix_mu=[0,0], fix_sig=[0,0], fix_a=False):
    """Fit the sample s with 2-gaussian mixture (histogram fit).
    Uses scipy.optimize.leastsq function. Parameters can be fixed but
    cannot be constrained in an interval.
    `p0`: initial guess or parameters
    `bins`: bins passed to `histogram()`
    `weights` optional weights for `histogram()`
    `fix_a`: if True fixes `a` to p0[4]
    `fix_mu`: tuple of bools. Whether to fix the means of the gaussians
    `fix_sig`: tuple of bools. Whether to fix the sigmas of the gaussians
    """
    assert size(p0) == 5
    fix = array([fix_mu[0], fix_sig[0], fix_mu[1], fix_sig[1], fix_a], 
                dtype=bool)
    p0 = array(p0)
    p0_free = p0[-fix]
    p0_fix = p0[fix]

    H = histogram(s, bins=bins, weights=weights, density=True)
    x, y = 0.5*(H[1][:-1] + H[1][1:]), H[0]
    assert x.size == y.size
    
    ## Fitting
    def err_func(p, x, y, fix, p_fix, p_complete): 
        p_complete[-fix] = p
        p_complete[fix] = p_fix
        return y - two_gauss_mix_pdf(x, p_complete)

    p_complete = zeros(5)
    p, v = leastsq(err_func, x0=p0_free, args=(x,y,fix,p0_fix,p_complete))
    p_new = zeros(5)
    p_new[-fix] = p
    p_new[fix] = p0_fix
    return reorder_parameters(p_new)
    
def two_gaussian_fit_hist_min(s, bounds=None, method='L-BFGS-B', 
        bins=r_[-0.5:1.5:0.001], weights=None,  p0=[0.2,1,0.8,1,0.3], 
        fix_mu=[0,0], fix_sig=[0,0], fix_a=False, verbose=False):
    """Fit the sample `s` with 2-gaussian mixture (histogram fit). [Bounded]
    Uses scipy.optimize.minimize allowing constrained minimization.
    `method`: one of the methods accepted by minimize()
    `bounds`: 5-element list of (min,max) values to constrain each of the 5
            parameters (can be used only with L-BFGS-B, TNC or SLSQP methods)
            If bounds are used, parameters cannot be fixed
    `p0`: initial guess or parameters
    `bins`: bins passed to `histogram()`
    `weights` optional weights for `histogram()`
    `fix_a`: if True fixes `a` to p0[4]
    `fix_mu`: tuple of bools. Whether to fix the means of the gaussians
    `fix_sig`: tuple of bools. Whether to fix the sigmas of the gaussians
    """
    assert size(p0) == 5
    fix = array([fix_mu[0], fix_sig[0], fix_mu[1], fix_sig[1], fix_a], 
                dtype=bool)
    p0 = array(p0)
    p0_free = p0[-fix]
    p0_fix = p0[fix]

    H = histogram(s, bins=bins, weights=weights, density=True)
    x, y = 0.5*(H[1][:-1] + H[1][1:]), H[0]
    assert x.size == y.size
    
    ## Fitting
    def err_func(p, x, y, fix, p_fix, p_complete): 
        p_complete[-fix] = p
        p_complete[fix] = p_fix
        return ((y - two_gauss_mix_pdf(x, p_complete))**2).sum()
    
    p_complete = zeros(5)
    res = minimize(err_func, x0=p0_free, args=(x,y,fix,p0_fix,p_complete),
                    method=method, bounds=bounds)
    if verbose: print(res)
    p_new = zeros(5)
    p_new[-fix] = res.x
    p_new[fix] = p0_fix
    return reorder_parameters(p_new)

def two_gaussian_fit_cdf(s, p0=[0.,.05,.6,.1,.5], fix_mu=[0,0], fix_sig=[0,0]):
    """Fit the sample s with two gaussians. Parameters can be fixed.
    """
    assert size(p0) == 5
    fix = array([fix_mu[0], fix_sig[0], fix_mu[1], fix_sig[1], 0], dtype=bool)
    p0 = array(p0)
    p0_free = p0[-fix]
    p0_fix = p0[fix]

    ## Empirical CDF
    ecdf = [sort(s), arange(0.5,s.size+0.5)*1./s.size]
    x, y = ecdf
    
    ## Analytical gaussian CDF
    gauss_cdf = lambda x, mu, sigma: 0.5*(1+erf((x-mu)/(sqrt(2)*sigma)))
    def two_gauss_mix_cdf(x, p):
        return p[4]*gauss_cdf(x,p[0],p[1])+(1-p[4])*gauss_cdf(x,p[2],p[3])

    ## Fitting the empirical CDF
    def err_func(p, x, y, fix, p_fix, p_complete): 
        p_complete[-fix] = p
        p_complete[fix] = p_fix
        return y - two_gauss_mix_cdf(x, p_complete)

    p_complete = zeros(5)
    p, v = leastsq(err_func, x0=p0_free, args=(x,y,fix,p0_fix,p_complete))
    p_new = zeros(5)
    p_new[-fix] = p
    p_new[fix] = p0_fix
    return reorder_parameters(p_new)

def test_two_gauss():
    m01 = 0.
    m02 = 0.6
    s01 = 0.05
    s02 = 0.1
    alpha = 0.
    p_real = [m01, s01, m02, s02, alpha]

    N = 500
    si1 = round(alpha*N)
    si2 = round((1-alpha)*N)
    s1 = R.normal(size=si1, loc=m01, scale=s01)
    s2 = R.normal(size=si2, loc=m02, scale=s02)
    s = r_[s1,s2]
    
    pc = two_gaussian_fit_cdf(s, fix_mu=[1,0], p0=[-0.01,0.05,0.5,0.2,0.4])
    ph = two_gaussian_fit_hist(s, fix_mu=[1,0], p0=[-0.01,0.05,0.5,0.2,0.4])
    pe = two_gaussian_fit_EM(s, fix_mu=[1,0], p0=[-0.01,0.05,0.5,0.2,0.4])
    
    hist(s, bins=40, normed=True)
    
    x = r_[s.min()-1:s.max()+1:200j]
    plot(x, a*normpdf(x,mu1,sig1), lw=2)
    plot(x, (1-a)*normpdf(x,mu2,sig2), lw=2)
    plot(x, two_gauss_mix_pdf(x, p0), lw=2)

    axvline(m01, lw=2, color='k', alpha=0.3)
    axvline(m02, lw=2, color='gray', alpha=0.3)
    axvline(mu1, lw=2, ls='--', color='k', alpha=0.3)
    axvline(mu2, lw=2, ls='--', color='gray', alpha=0.3)
    axvline(mu1h, lw=2, ls='--', color='r', alpha=0.3)
    axvline(mu2h, lw=2, ls='--', color='r', alpha=0.3)

def compare_two_gauss():
    m01 = 0.
    m02 = 0.5
    s01 = 0.08
    s02 = 0.15
    alpha = 0.7
    p_real = [m01, s01, m02, s02, alpha]

    N = 1000
    si1 = round(N*alpha)
    si2 = round((1-alpha)*N)
    
    p0 = [-0.01,0.05,0.6,0.2,0.4]
    fix_mu = [0,0]

    n = 500
    PC, PH, PE = zeros((n,5)), zeros((n,5)), zeros((n,5))
    for i in xrange(n):
        s1 = R.normal(size=si1, loc=m01, scale=s01)
        s2 = R.normal(size=si2, loc=m02, scale=s02)
        s = r_[s1,s2]
        
        pc = two_gaussian_fit_cdf(s, fix_mu=fix_mu, p0=p0)
        ph = two_gaussian_fit_hist(s, fix_mu=fix_mu, p0=p0)
        pe = two_gaussian_fit_EM(s, fix_mu=fix_mu, p0=p0)
        
        PC[i], PH[i], PE[i] = pc, ph, pe
    
    Label = ['Mu1', 'Sig1', 'Mu2', 'Sig2', 'Alpha']
    ftype = 'png'
    
    for i in range(5):
        figure()
        title(Label[i])
        vmin = min([PC[:,i].min(), PH[:,i].min(), PE[:,i].min()])
        vmax = max([PC[:,i].max(), PH[:,i].max(), PE[:,i].max()])
        b = r_[vmin:vmax:80j]
        if vmax == vmin: b = r_[vmin-.1:vmax+.1:200j]
        hist(PC[:,i], bins=b, alpha=0.3, label='CDF')
        hist(PH[:,i], bins=b, alpha=0.3, label='Hist')
        hist(PE[:,i], bins=b, alpha=0.3, label='EM')
        legend(loc='best')
        axvline(p_real[i], color='k', lw=2)
        #savefig('Two-gaussian Fit Comp - %s.png' % Label[i])

def gaussian2d_fit(sx, sy, guess=[0.5,1]):
    """2D-Gaussian fit of samples S using a fit to the empirical CDF."""
    assert sx.size == sy.size
    
    ## Empirical CDF
    ecdfx = [sort(sx), arange(0.5,sx.size+0.5)*1./sx.size]
    ecdfy = [sort(sy), arange(0.5,sy.size+0.5)*1./sy.size]
    
    ## Analytical gaussian CDF
    gauss_cdf = lambda x, mu, sigma: 0.5*(1+erf((x-mu)/(sqrt(2)*sigma)))
    
    ## Fitting the empirical CDF
    fitfunc = lambda p, x: gauss_cdf(x, p[0], p[1])
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    px,v = leastsq(errfunc, x0=guess, args=(ecdfx[0],ecdfx[1]))
    py,v = leastsq(errfunc, x0=guess, args=(ecdfy[0],ecdfy[1]))
    print "2D Gaussian CDF fit", px, py
    
    mux, sigmax = px[0], px[1]
    muy, sigmay = py[0], py[1]
    return mux, sigmax, muy, sigmay

def test_gaussian2d_fit():
    mx0 = 0.1
    my0 = 0.9
    sigx0 = 0.4
    sigy0 = 0.25

    Size = 500
    sx = R.normal(size=Size, loc=mx0, scale=sigx0)
    sy = R.normal(size=Size, loc=my0, scale=sigy0)

    mux, sigmax, muy, sigmay = gaussian2d_fit(sx, sy)

    plot(sx, sy, 'o', alpha=0.2, mew=0)

    X,Y = mgrid[sx.min()-1:sx.max()+1:200j, sy.min()-1:sy.max()+1:200j]
    
    def gauss2d(X,Y, mx, my, sigx, sigy):
        return exp(-((X-mx)**2)/(2*sigx**2))*exp(-((Y-my)**2)/(2*sigy**2))
    
    contour(X,Y,gauss2d(X,Y,mux,muy,sigmax,sigmay))

    plot(mx0,my0, 'ok', mew=0, ms=10)
    plot(mux,muy, 'x', mew=2, ms=10, color='green')


def two_gaussian2d_fit(sx, sy, guess=[0.5,1]):
    """2D-Gaussian fit of samples S using a fit to the empirical CDF."""
    ## UNFINISHED (I have 2 alphas usign the xy projections)
    assert sx.size == sy.size
    
    ## Empirical CDF
    ecdfx = [sort(sx), arange(0.5,sx.size+0.5)*1./sx.size]
    ecdfy = [sort(sy), arange(0.5,sy.size+0.5)*1./sy.size]
    
    ## Analytical gaussian CDF
    gauss_cdf = lambda x, mu, sigma: 0.5*(1+erf((x-mu)/(sqrt(2)*sigma)))

    gauss2d_cdf = lambda X,Y,mx,sx,my,sy: gauss_cdf(X,mx,sx)*gauss_cdf(Y,my,sy)

    two_cdf = lambda x, m1, s1, m2, s2, a:\
        a*gauss_cdf(x,m1,s1)+(1-a)*gauss_cdf(x,m2,s2)

    two2d_cdf = lambda X,Y, mx1, sx1, mx2, sx2, my1, sy1, my2, sy2, a:\
        a*gauss2d_cdf(X,Y,mx1,sx1,my1,sy1)+(1-a)*gauss_cdf(X,Y,mx2,sx2,my2,sy2)

    ## Fitting the empirical CDF
    fitfunc = lambda p, x: two_cdf(x, *p)
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    fitfunc2d = lambda p, X,Y: two2d_cdf(X,Y, *p)
    errfunc2d = lambda p, X,Y,Z: fitfunc2d(p, X,Y) - Z

    px,v = leastsq(errfunc, x0=guess, args=(ecdfx[0],ecdfx[1]))
    py,v = leastsq(errfunc, x0=guess, args=(ecdfy[0],ecdfy[1]))
    print "2D Two-Gaussians CDF fit", px, py
    
    mux1, sigmax1, mux2, sigmax2, alphax = px
    muy1, sigmay1, muy2, sigmay2, alphay = py
    return mu1, sigma1, mu2, sigma2, alpha

def test_gaussian_fit():
    m0 = 0.1
    s0 = 0.4
    size = 500

    s = R.normal(size=size, loc=m0, scale=s0)
    #s = s[s<0.4]
    mu, sig = gaussian_fit(s)
    mu1, sig1 = S.norm.fit(s)
    mu2, sig2 = gaussian_fit_ml(s)
    
    print "ECDF ", mu, sig
    print "ML         ", mu1, sig1
    print "ML (manual)", mu2, sig2
    
    H = histogram(s, bins=20, density=True)
    h = H[0]
    bw = H[1][1] - H[1][0]
    #bins_c = H[1][:-1]+0.5*bw
    bar(H[1][:-1], H[0], bw, alpha=0.3)

    x = r_[s.min()-1:s.max()+1:200j]
    plot(x, normpdf(x,m0,s0), lw=2, color='grey')    
    plot(x, normpdf(x,mu,sig), lw=2, color='r', alpha=0.5)
    plot(x, normpdf(x,mu1,sig1), lw=2, color='b', alpha=0.5)
    
if __name__ == '__main__':
    #compare_two_gauss()
    #test_gaussian2d_fit()
    #test_gaussian_fit()
    #show()
    pass

