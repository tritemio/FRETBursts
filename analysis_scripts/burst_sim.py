"""
Simulate a BSD and try different fit functions.

NOTE: Old version, the notebooks contain much more comprehensive simulations
      and fit comparisons.
"""


import scipy.stats as SS
import bt_fit

## Total number of burts to draw
N_bursts = 1e3

## Tau of the exponential distribution of burst sizes
burst_size_scale = 20

## Slope of the width vs size relation, and width in case size zero
width_size_slope = 10e-6
width_size0 = 0.5e-3

## Ratio std(width)/mean(width)
width_std_mean_ratio = 0.2

## BG rates
bg_d_rate = 9e3
bg_a_rate = 7e3

## Bleedthrough
BT = 0.05

#  - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Draw 1: Donor burst size
Sd = SS.expon(loc=0, scale=burst_size_scale).rvs(N_bursts)
#Sd = floor(Sd)

## Draw 2: Burst widths
mean_width = Sd*width_size_slope + width_size0
std_width = width_std_mean_ratio*mean_width
w = r_[[SS.norm(loc=mean_width[i], scale=std_width[i]).rvs() 
    for i in xrange(size(Sd))]]

## Draw 3: BG in donor and accept channels
bg_d_mean = bg_d_rate*w
bg_a_mean = bg_a_rate*w
bg_d_rv = r_[[SS.poisson(bgi).rvs() for bgi in bg_d_mean]] 
bg_a_rv = r_[[SS.poisson(bgi).rvs() for bgi in bg_a_mean]] 

## Draw 4: Bleedthrough in the A channel
bt_n = zeros(N_bursts)
bt_n[Sd>0] = r_[[SS.poisson(BT*x).rvs() for x in Sd[Sd>0]]]

## Build the burst data
nd = Sd+bg_d_rv
na = bg_a_rv+bt_n
bg_d, bg_a = bg_d_mean, bg_a_mean

## BT fitting
bt_lr = bt_fit.fit_LR_s(nd,na,bg_d,bg_a)
bt_ls = bt_fit.fit_LS_s(nd,na,bg_d,bg_a)
bt_ml = bt_fit.fit_ML_s(nd,na,bg_d,bg_a)

print "Max-likelihood: %.4f%%, error %.3f%%" %(bt_ml*100, (BT-bt_ml)/BT*100)
print "Linear regress: %.4f%%, error %.3f%%" %(bt_ls*100, (BT-bt_ls)/BT*100)
