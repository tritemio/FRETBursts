"""
## Only intensity profile
Cy3B_ND3_accept_xm_ym_0.dat
Cy3B_ND4_donor_xm_ym_1.dat

## BT measurements at high concentration
12d_DO_ND34__xm_ym_all_2.dat
12d_DO_ND34__xm_ym_DCR_3.dat

12d_DO_ND34__xm_ym_1off_4.dat
12d_DO_ND34__xm_ym_2off_5.dat
12d_DO_ND34__xm_ym_3off_6.dat
12d_DO_ND34__xm_ym_4off_7.dat
12d_DO_ND34__xm_ym_5off_8.dat
12d_DO_ND34__xm_ym_6off_9.dat
12d_DO_ND34__xm_ym_7off_10.dat
12d_DO_ND34__xm_ym_8off_11.dat

12d_DO_ND34__xm_ym_all_12.dat
12d_DO_ND34__xm_ym_all_steer_13.dat

12d_DO_ND34__xm_ym_1off_steer_14.dat
12d_DO_ND34__xm_ym_2off_steer_15.dat
12d_DO_ND34__xm_ym_3off_steer_16.dat
12d_DO_ND34__xm_ym_4off_steer_17.dat
12d_DO_ND34__xm_ym_5off_steer_18.dat
12d_DO_ND34__xm_ym_6off_steer_19.dat
12d_DO_ND34__xm_ym_7off_steer_20.dat
12d_DO_ND34__xm_ym_8off_steer_21.dat

12d_DO_ND34__xm_ym_all_steer_22.dat
12d_DO_ND34__xm_ym_all_23.dat

12d_DO_ND34__xm_ym_all_24.dat
12d_DO_ND34__xm_ym_all_25.dat

## --> From here single molecule data:
12d_DO_420p_320mW_steering_26.dat
12d_DO_420p_320mW_27.dat

12d_DO_150p_320mW_28.dat
12d_DO_150p_320mW_steering_29.dat
12d_DO_150p_320mW_steering_35_4-2_60_30.dat

12d_DO_150pNo2_320mW_31.dat
12d_DO_150pNo2_320mW_steering_32.dat
"""

f1 = '12d_DO_ND34__xm_ym_1off_4.dat'
f2 = '12d_DO_ND34__xm_ym_2off_5.dat'
f3 = '12d_DO_ND34__xm_ym_3off_6.dat'
f4 = '12d_DO_ND34__xm_ym_4off_7.dat'
f5 = '12d_DO_ND34__xm_ym_5off_8.dat'
f6 = '12d_DO_ND34__xm_ym_6off_9.dat'
f7 = '12d_DO_ND34__xm_ym_7off_10.dat'
f8 = '12d_DO_ND34__xm_ym_8off_11.dat'
F = [f1,f2,f3,f4,f5,f6,f7,f8]

fs1 = '12d_DO_ND34__xm_ym_1off_steer_14.dat'
fs2 = '12d_DO_ND34__xm_ym_2off_steer_15.dat'
fs3 = '12d_DO_ND34__xm_ym_3off_steer_16.dat'
fs4 = '12d_DO_ND34__xm_ym_4off_steer_17.dat'
fs5 = '12d_DO_ND34__xm_ym_5off_steer_18.dat'
fs6 = '12d_DO_ND34__xm_ym_6off_steer_19.dat'
fs7 = '12d_DO_ND34__xm_ym_7off_steer_20.dat'
fs8 = '12d_DO_ND34__xm_ym_8off_steer_21.dat'
FS = [fs1,fs2,fs3,fs4,fs5,fs6,fs7,fs8]

fa1 = '12d_DO_ND34__xm_ym_all_2.dat'
f_dcr = '12d_DO_ND34__xm_ym_DCR_3.dat'

fa2 = '12d_DO_ND34__xm_ym_all_12.dat'
fsa1 = '12d_DO_ND34__xm_ym_all_steer_13.dat'

fsa2 = '12d_DO_ND34__xm_ym_all_steer_22.dat'
fa3 = '12d_DO_ND34__xm_ym_all_23.dat'

def rate_ph_total(fname):
    dir_ = 'C:/Data/Antonio/data/2013-05-09/'
    d = Data(fname=dir_+fname, clk_p=clk_p, nch=8, BT=0., gamma=1.)
    d.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
    rates = array([1.*ph.size/(ph[-1]-ph[0]) for ph in d.ph_times_m])
    rates /= d.clk_p
    return rates

def rate_ph_da(fname):
    dir_ = 'C:/Data/Antonio/data/2013-05-09/'
    d = Data(fname=dir_+fname, clk_p=12.5e-9, nch=8, BT=0., gamma=1.)
    d.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
    rates_a = array([1.*ph[a].size/(ph[a][-1]-ph[a][0])
            for ph,a in zip(d.ph_times_m,d.A_em)])
    rates_d = array([1.*ph[-a].size/(ph[-a][-1]-ph[-a][0])
            for ph,a in zip(d.ph_times_m,d.A_em)])
    rates_d /= d.clk_p
    rates_a /= d.clk_p
    return rates_d, rates_a

#R = [rate_ph_total(f) for f in F]
#RS = [rate_ph_total(f) for f in FS]

#Rda = array([rate_ph_da(f) for f in F])
RSda = array([rate_ph_da(f) for f in FS])

rt = rate_ph_da(fa1)
r_dcr = rate_ph_da(f_dcr)


def bt_from_rates(R, rt):
    """R must contain D and A rates for 8 measurents, each with a spot off.
    rt contains the D and A rates for when all spots are ON
    """
    assert R.shape == (8L,2L,8L)
    BT = []
    for i,s in enumerate(R):
        don = rt[0][i] - s[0][i]
        acc = rt[1][i] - s[1][i]
        BT.append(acc/don)
    return BT

BT1 = bt_from_rates(Rda)
BT2 = bt_from_rates(RSda) # meas. using "beam steering"
