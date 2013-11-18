import scipy.optimize as O
import scipy.stats as S
import numpy.random as R
from scipy.special import erf, erfc
from scipy.optimize import leastsq

def exp_gauss_pdf(x, mu, sig, lamb):
    return 0.5*lamb*exp(0.5*lamb*(2*mu+lamb*sig**2-2*x)) * \
            erfc((mu+lamb*sig**2-x)/(sqrt(2)*sig))
