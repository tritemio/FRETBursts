"""
Simple test for burst search and burst fuse functions.
"""

import pickle
ip = get_ipython()
ip.magic("run -i ../burstsearch.py")

ph = pickle.load(open("ph.pickle"))

clk_p = 12.5e-9
ph_ms = ph*clk_p*1e3

h = histogram(ph, bins=r_[:7000:1]/(12.5e-9*1e3))
plot(h[1][:-1],h[0])

T=round(0.7/clk_p/1e3)
mb = ba(ph, L=35, m=10, T=T)

for b in mb:
    axvline(b[itstart], color='k')
    axvline(b[itend], color='r')
    axvspan(b[itstart], b[itend], color='k', alpha=0.2)

# Burst 21 and 22 are superimposed

mbf = b_fuse(ph, mb, ms=0, clk_p=clk_p)

figure()
plot(h[1][:-1],h[0])
for b in mbf:
    axvline(b[itstart], color='k')
    axvline(b[itend], color='r')
    axvspan(b[itstart], b[itend], color='k', alpha=0.2)


