# encoding: utf-8
#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
WARNING: Plot function ATTIC! Functions here are broken!

Here there are function originally placed in burst_plot.py that became
broken and I didn't had the time or the need to update. They live here until
I decide to fix or delete them.
"""

##
#  Other plot wrapper functions
#

def wplot(*args, **kwargs):
    AX, s = dplot_8ch(*args, **kwargs)
    kwargs.update(AX=AX)
    q = gs.mToolQT(gcf(), dplot_8ch, *args, **kwargs)
    return AX, q

def splot(d, fun=scatter_width_size,
        scroll=False, pgrid=True, figsize=(10, 8), nosuptitle=False, ax=None,
        scale=True, **kwargs):
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    for i in xrange(d.nch):
        try:
            if i == 0 and not nosuptitle: fig.set_title(d.status())
        except:
            print("WARNING: No title in plots.")
        ax.grid(pgrid)
        fun(d, i, **kwargs)
    s = None
    if scroll: s = ScrollingToolQT(fig)
    return ax, s

##
#  Other misc plot functions
#
def bplot(d, ich, b_index,  ph0=True, pad=0):
    """Plot photons in a burst as vertical bars. Burst: d.mburst[ich][b_index].
    """
    br = bl.b_irange(d.mburst[ich], b_index, pad=pad)
    accept = d.A_em[ich][br]
    donor = -accept
    ph = d.ph_times_m[ich][br]
    if ph0: ph -= ph[0]
    dtime = (ph[donor])*d.clk_p*1e6
    atime = (ph[accept])*d.clk_p*1e6
    plt.vlines(dtime,0,1, lw=2, color='g', alpha=0.8)
    plt.vlines(atime,0,1, lw=2, color='r', alpha=0.8)
    #plot(dtime, ones_like(ph[donor]), '^', color='g', alpha=0.5)
    #plot(atime, -ones_like(ph[accept]), 'o', color='r', alpha=0.5)
    xlabel("Time (us)")
    nd, na, nt = donor.sum(), accept.sum(), ph.size
    E = float(na)/(nt)
    title("#ph = %d, #D-ph = %d, #A-ph = %d, E = %.2f" % (nt,nd,na,E))
    plt.ylim(-10,10)


def bg_legend_8ch(d):
    ax = gca()
    L = ax.get_lines()[1::2]
    for i,l in enumerate(L):
        ich = i/3
        x = i%3
        s = ['Tot', 'D', 'A']
        r = [d.rate_m[ich], d.rate_dd[ich], d.rate_ad[ich]]
        l.set_label("CH%d, %s %d cps" % (ich+1, s[x], r[x]))
    ax.legend()
    plt.draw()

def bg_legend_alex(d):
    ax = gca()
    L = ax.get_lines()[1::2]
    for i,l in enumerate(L):
        ich = i/4
        x = i%4
        s = ['Tot', 'DD', 'AD', 'AA']
        r = [d.rate_m[ich], d.rate_dd[ich], d.rate_ad[ich], d.rate_aa[ich]]
        l.set_label("CH%d, %s %d cps" % (ich+1, s[x], r[x]))
    ax.legend()
    plt.draw()
