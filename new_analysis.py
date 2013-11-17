d7.update_bt(BT=0)
d7.burst_search_t(L=10, m=10, P=None, F=6, dither=dither, ph_sel='DA')
ds = Sel(d7, select_bursts_nda, th1=20, gamma1=1.)
ds.fit_E_generic(E2=0.2, fit_fun=gaussian_fit_hist, weights='size')*100
E0 = ds.E_fit.copy()
chi0 = calc_chi(E0, E0.mean())

d7.burst_search_t(L=10, m=10, P=None, F=6, dither=dither, ph_sel='A')
ds = Sel(d7, select_bursts_nda, th1=20, gamma1=.45)
ds.fit_E_generic(E1=0.5, fit_fun=gaussian_fit_hist, weights='size')*100
E1 = ds.E_fit.copy()
chi1 = calc_chi(E1, E1.mean())

chi_m = vstack([chi0,chi1]).mean(axis=0)
chi_m
plot(chi0)
plot(chi1)
plot(chi_m)

# BT method1
E0_m = E0.mean()
k_bt = E0_m/(1-E0_m)

# BT method 2
k_bt = mean(E0/(1-E0))

k_bt

d7.update_bt(BT=k_bt)
d7.update_chi_ch(chi_ch=chi_m)
d7.burst_search_t(L=10, m=10, P=None, F=6, dither=dither, ph_sel='DA')

def delta(x): return x.max()-x.min()

d7.burst_search_t(L=10, m=10, P=None, F=6, dither=dither, ph_sel='D')
ds = Sel(d7, select_bursts_nda, th1=20, gamma1=1.)
ds.fit_E_generic(E2=0.2, fit_fun=gaussian_fit_hist, weights='size')*100
ds.fit_E_generic(E2=0.2, fit_fun=gaussian_fit_hist, weights='size')
E0_final = ds.E_fit.copy() 
E0_final*100
delta(E0_final)*100

d7.burst_search_t(L=10, m=10, P=None, F=6, dither=dither, ph_sel='DA')
ds = Sel(d7, select_bursts_nda, th1=20, gamma1=.45)
ds.fit_E_generic(E1=0.5, fit_fun=gaussian_fit_hist, weights='size')*100
E1_final = ds.E_fit.copy()
E1_final*100
delta(E1_final)*100


correct_E(E1, gamma)  # don't apply chi_m, it's already in the fit

ds.BT
ds.gamma
ds.chi_ch
