# -*- coding: utf-8 -*-
"""
Test the effect of beam-steering on the spot-pattern quality.
"""

CH = 0 # 0 for donor, 1 for accept
sel = ["donor", "accept"]
meas = "_R2"

fa = "2013-05-15/Cy3B_%s%s_all_?.dat" % (sel[CH], meas)
fsa = "2013-05-15/Cy3B_%s%s_all_steer_?.dat" % (sel[CH], meas)

fc = "2013-05-15/Cy3B_%s%s_center_off_?.dat" % (sel[CH], meas)
fsc = "2013-05-15/Cy3B_%s%s_center_off_steer_?.dat" % (sel[CH], meas)

fsides = "2013-05-15/Cy3B_%s%s_sides_off_?.dat" % (sel[CH], meas)
fssides = "2013-05-15/Cy3B_%s%s_sides_off_steer_?.dat" % (sel[CH], meas)

f_dcr = "2013-05-15/Cy3B_%s_dcr_?.dat" % sel[CH]

from glob import glob
try:
    ip = get_ipython()
except:
    ip = _ip
ip.magic("run -i burst_selection.py")


def rate_ph_da(fname, kcps=True):
    fn = glob(data_dir+fname+'*')[0]
    
    d = Data(fname=fn, clk_p=12.5e-9, nch=8, BT=0., gamma=1.)
    d.load_multispot(bytes_to_read=-1, swap_D_A=True)
    rates_a, rates_d = zeros(8), zeros(8)
    for i,(ph,a) in enumerate(zip(d.ph_times_m,d.A_em)):
        if a.any():
            rates_a[i] = 1.*ph[a].size/(ph[a][-1]-ph[a][0])
    for i,(ph,a) in enumerate(zip(d.ph_times_m,d.A_em)):
        if (-a).any(): 
            rates_d[i] = 1.*ph[-a].size/(ph[-a][-1]-ph[-a][0])
    rates_d /= d.clk_p
    rates_a /= d.clk_p
    if kcps:
        rates_d /= 1e3
        rates_a /= 1e3
    return rates_d, rates_a

## Compute donor rates
if 1:
    RA = rate_ph_da(fa)
    RSA = rate_ph_da(fsa)
    
    RC = rate_ph_da(fc)
    RSC = rate_ph_da(fsc)
    
    RSides = rate_ph_da(fsides)
    RSSides = rate_ph_da(fssides)

    R_dcr = rate_ph_da(f_dcr)
    
    RA,RSA,RC,RSC,RSides,RSSides,R_dcr = \
        RA[CH],RSA[CH],RC[CH],RSC[CH],RSides[CH],RSSides[CH],R_dcr[CH]
    
    for R in [RA,RSA,RC,RSC,RSides,RSSides]:
        R -= R_dcr

max_cps, max_ch = RA.max(), RA.argmax()+1
max_cps_s, max_ch_s = RSA.max(), RSA.argmax()+1
min_cps, min_ch = RA.min(), RA.argmin()+1
min_cps_s, min_ch_s = RSA.min(), RSA.argmin()+1
Ratio_max_min = max_cps/min_cps
Ratio_max_min_s = max_cps_s/min_cps_s

center_cps, center_cps_off = RA[4], RC[4]
center_cps_s, center_cps_off_s = RSA[4], RSC[4]
Ratio_center = center_cps/center_cps_off
Ratio_center_s = center_cps_s/center_cps_off_s

side1_cps, side1_cps_off = RA[7], RSides[7]
side1_cps_s, side1_cps_off_s = RSA[7], RSSides[7]
Ratio_side1 = side1_cps/side1_cps_off 
Ratio_side1_s = side1_cps_s/side1_cps_off_s

side2_cps, side2_cps_off = RA[0], RSides[0]
side2_cps_s, side2_cps_off_s = RSA[0], RSSides[0]
Ratio_side2 = side2_cps/side2_cps_off 
Ratio_side2_s = side1_cps_s/side1_cps_off_s

print "\n => CH: %s MEASUREMENT %s" % (sel[CH], meas)
print "Max: %3d (CH%d)   Min: %3d (CH%d)   Ratio: %4.2f" %\
        (max_cps,max_ch,min_cps,min_ch,Ratio_max_min)
print "Max: %3d (CH%d)   Min: %3d (CH%d)   Ratio: %4.2f" %\
        (max_cps_s,max_ch_s,min_cps_s,min_ch_s,Ratio_max_min_s)

print "            ON  OFF  Ratio"
print "Center:    %3d  %3d   %4.2f" %(center_cps, center_cps_off, Ratio_center)
print "Center(s): %3d  %3d   %4.2f" %(center_cps_s, center_cps_off_s, Ratio_center_s)

print "Side1:     %3d  %3d   %4.2f" %(side1_cps, side1_cps_off, Ratio_side1)
print "Side1(s):  %3d  %3d   %4.2f" %(side1_cps_s, side1_cps_off_s, Ratio_side1_s)

print "Side2:     %3d  %3d   %4.2f" %(side2_cps, side2_cps_off, Ratio_side2)
print "Side2(s):  %3d  %3d   %4.2f" %(side2_cps_s, side2_cps_off_s, Ratio_side2_s)

plot(RA); plot(RSA); ylim(0)
