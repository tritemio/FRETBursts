"""
Compute Mean burst size (BS), Number of bursts (NB) and number of photon (PH)
for different values of F.

"""

COLORS = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange']
CH = r_[1:9]

m = 5
FF = r_[2:10.1:0.5]
BS, NB, PH = [], [], []
for F in FF:
    d.burst_search_t(L=m,m=m,F=F,P=None, nofret=0)
    #d.cal_ph_num()
    BS.append([mean(nt) for nt in d.nt])
    NB.append(d.num_bu())
    PH.append([nt.sum() for nt in d.nt])
BS, NB, PH = array(BS), array(NB), array(PH)


# Plot mean burst size vs F
figure()
title("Mean BS vs F (m=%d)"%m)
plot(FF, BS, lw=2)
legend(['CH%d' % ch for ch in CH], frameon=0, ncol=3, loc='best')
xlabel('F'); grid(True)

# Mark the maximum posision with a vertical line
for c, ch in zip(COLORS,CH):
    ax = axvline(FF[BS[:,ch-1].argmax()], ymin=(1-ch/8.), lw=10-ch, alpha=0.5, color=c)


savefig("temp/12d_BS_vs_F_m%02d.png" % m, bbox_inches='tight')
