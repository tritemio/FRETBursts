from time import sleep

class AlexHist:
    def __init__(self, d, bin_step=0.02):
        self.d = d
        self.bin_step = bin_step
        self.compute()
    def compute(self):
        H, E_bins, S_bins = ES_histo(self.d.E[0],self.d.S[0], self.bin_step)
        self.H, self.E_bins, self.S_bins = H, E_bins, S_bins
        self.E_ax = self.E_bins[:-1] +0.5*self.bin_step
        self.S_ax = self.S_bins[:-1] +0.5*self.bin_step

def ES_histo(E,S, bin_step=0.05, E_bins=None, S_bins=None):
    if E_bins is None: E_bins = arange(-0.2,1.2+1e-4,bin_step)
    if S_bins is None: S_bins = arange(-0.2,1.2+1e-4,bin_step)
    H, E_bins, S_bins = histogram2d(E,S,bins=[E_bins,S_bins])
    return H, E_bins, S_bins

def alexhist(d, vmin=2, vmax=0, interp='bicubic', 
        cmap='hot', under_color='white', over_color='white', ax=None,
        scatter=True, scatter_ms=3, scatter_color='orange', scatter_alpha=0.2):
    
    AH, E_bins,S_bins, E_ax,S_ax = d.AH[i], d.E_bins,d.S_bins, d.E_ax,d.S_ax

    if ax is None: ax = subplot(111)
    colormap = get_cmap(cmap)
    if vmax < vmin:
        E_range = (E_bins > 0.4)*(E_bins < 0.8)
        vmax = AH[E_range,:].max()
    if scatter:
        ax.plot(d.E[0],d.S[0], 'o', mew=0, ms=scatter_ms, alpha=scatter_alpha, 
                color=scatter_color)
    im = ax.imshow(AH[:,::-1].T, interpolation=interp,
            extent=(E_bins[0],E_bins[-1],S_bins[0],S_bins[-1]), 
            vmin=vmin, vmax=vmax, cmap=colormap)
    im.cmap.set_under(under_color)
    im.cmap.set_over(over_color)
    ax.figure.colorbar(im)
    ax.set_xlim(-0.2,1.2); ax.set_ylim(-0.2,1.2)
    ax.set_xlabel('E'); ax.set_ylabel('S')
    return ax

def alexhist_imshow(alex_hist, vmin=2, vmax=0, interp='bicubic', 
        cmap='hot', under_color='white', over_color='white', ax=None,
        scatter=True, scatter_ms=3, scatter_color='orange', scatter_alpha=0.2):
    H, E_bins, S_bins, E_ax, S_ax = alex_hist.H, alex_hist.E_bins, \
            alex_hist.S_bins, alex_hist.E_ax, alex_hist.S_ax

    E, S = alex_hist.d.E[0], alex_hist.d.S[0]

    if ax is None: ax = subplot(111)
    colormap = get_cmap(cmap)
    if vmax < vmin:
        E_range = (E_bins > 0.4)*(E_bins < 0.8)
        vmax = H[E_range,:].max()
    if scatter:
        ax.plot(d.E[0],d.S[0], 'o', mew=0, ms=scatter_ms, alpha=scatter_alpha, 
                color=scatter_color)
    im = ax.imshow(H[:,::-1].T, interpolation=interp,
            extent=(E_bins[0],E_bins[-1],S_bins[0],S_bins[-1]), 
            vmin=vmin, vmax=vmax, cmap=colormap)
    im.cmap.set_under(under_color)
    im.cmap.set_over(over_color)
    ax.figure.colorbar(im)
    ax.set_xlim(-0.2,1.2); ax.set_ylim(-0.2,1.2)
    ax.set_xlabel('E'); ax.set_ylabel('S')
    return ax

def alexhist_contour(alex_hist, vmin=1, vmax=60, n_levels=8, 
        cmap='BuPu', under_color='white', over_color='white', ax=None):
    H, E_bins, S_bins, E_ax, S_ax = alex_hist.H, alex_hist.E_bins, \
            alex_hist.S_bins, alex_hist.E_ax, alex_hist.S_ax

    if ax is None: ax = subplot(111)
    levels = arange(vmin,vmax+1, (vmax-vmin)/n_levels)
    colormap = get_cmap(cmap)
    ct = ax.contour(E_ax, S_ax, H.T, levels=levels, alpha=1, 
            cmap=colormap, zorder=2, linewidths=(2,))
    ax.figure.colorbar(ct)
    ax.set_xlabel('E'); ax.set_ylabel('S')
    return ax

def alexhist_contourf(alex_hist, vmin=0, vmax=12, n_levels=8, 
        cmap='BuPu', under_color='white', over_color='white', ax=None):
    H, E_bins, S_bins, E_ax, S_ax = alex_hist.H, alex_hist.E_bins, \
            alex_hist.S_bins, alex_hist.E_ax, alex_hist.S_ax

    if ax is None: ax = subplot(111)
    levels = arange(vmin,vmax+1, (vmax-vmin)/n_levels)
    colormap = get_cmap(cmap)
    cf = ax.contourf(E_ax, S_ax, H.T, levels=levels, alpha=1, 
            cmap=colormap, zorder=2)
    ax.figure.colorbar(cf)
    ax.set_xlabel('E'); ax.set_ylabel('S')
    return ax

#AH = AlexHist(d)
#alexhist_imshow(AH, vmax=1000, interp='nearest')
