"""
Functions to simulate burst sizes.
"""

from scipy.stats import poisson
from scipy.misc import factorial
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar, leastsq
from scipy.special import gammaln

def E0_sim(Sd, BGa):
    """Simulation of D-only peak due to fixed donor values Sd, and bg in A BGa
    """
    assert size(Sd) == size(BGa)
    Sd, BGa = array(Sd, ndmin=1), array(BGa, ndmin=1)
    
    ## Continum expansion of a poisson PMF (shifted to have mean=0)
    poiss_pdf_m = lambda x,l: (exp(-l)*(l**(x+l)))/(factorial(x+l))
    
    ## Axis of possible BG-corrected A values (Poisson PMF domain)
    ## This should be big enough for all the values in BGa
    A_ax = r_[:BGa.max()+5*sqrt(BGa.max())+1:0.125*0.5] - BGa.max()
    
    ## Build a common E axis for intepolation
    E_ax = r_[-1:2:0.005]

    P_tot = zeros(E_ax.size)
    for sd,bga in zip(Sd, BGa):
        ## Possible E values
        E = A_ax/(A_ax+sd)

        ## Probability for each A_ax (and therefore E) value
        P = poiss_pdf_m(A_ax, bga)
        
        ## Reintepolate the E values on a common E axis
        P_i = interp1d(E, P, fill_value=0, bounds_error=False)(E_ax)
        
        ## Add current (sd,bga) value to the total histogram
        P_tot += P_i
    return E_ax, P_tot



