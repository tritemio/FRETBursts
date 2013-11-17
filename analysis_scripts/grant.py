def plot_mburstm_8ch_global(D, fun=hist_fret2, pgrid=True, figsize=(12,9), 
        **kwargs):
    fig, AAX = subplots(8,5,figsize=figsize, sharex=True, sharey=True)
    for AX,d in zip(AAX.T,D):
        for i, ax in enumerate(AX):
            b = d.mburst[i] if hasattr(d, 'mburst') else None   
            if (not fun in [timetrace, ratetrace]) and b.size is 0: continue
            sca(ax)
            fun(i,b,d,**kwargs)
            #text(0.4,0.7, "CH%d" % (i+1), transform = gca().transAxes,
            #    fontsize='large')
            ax.set_xlabel(''); ax.set_ylabel('')
            ax.set_xticks([0,0.5,1]); ax.set_yticks(r_[0:10])
            ax.tick_params(width=1.2, length=6)
            #ax.autoscale(enable=True,axis='y')
        setp([a.get_xticklabels() for a in AX[:-1]], visible=False)
        setp([a.get_yticklabels() for a in AX], visible=False)
    subplots_adjust(left=0.06, right=0.96, top=0.94, bottom=0.05)
    fig.subplots_adjust(hspace=0, wspace=0)
    setp([a.get_yticklabels() for a in AAX[:,0]], visible=True)
    for ax,sa in zip(AAX[0,:],['7bp','12bp','17bp','22bp','27bp']):
        ax.set_title(sa)
    for ich,ax in enumerate(AAX[:,0]): ax.set_ylabel("CH%d"%(ich+1))
    return AX

#AX = plot_mburstm_8ch_global([dfs7, dfs12, dfs17, dfs22,dfs27],bins=r_[-0.2:1.2:0.05]+0.025, normed=True)
#ylim(0,4.8)


if 0:
    rcParams['xtick.labelsize']='medium'
    rcParams['ytick.labelsize']='medium'
    AX, s = plot_mburstm_8ch(dfs12, hist_fret, bins=r_[-0.2:1.2:0.05]+0.025,
            pgrid=False)
    AX[-1].set_xticks([0,0.2,0.4,0.6,0.8,1])
    for ax in AX: ax.set_yticks(r_[0:7:2])
    AX[-1].set_xlabel('E', fontsize='large')
    tight_layout()
    fig = gcf()
    fig.subplots_adjust(hspace=0)
    for ax in AX: ax.tick_params(width=1.2, length=6) 
    #ax, s = plot_mburstm_8ch_super(dfs, hist_fret, bins=r_[-0.2:1.21:0.05])
    #plot_mburstm_8ch_global([dfs, dfs, dfs, dfs, dfs])

if 0:
    fontsize = 'large'
    rcParams['xtick.labelsize']=fontsize
    rcParams['ytick.labelsize']=fontsize
    AX,s = dplot(dfs12, hist_fret, figsize=(5, 5.975), 
                 bins=r_[-0.2:1.2:0.05]+0.025, pgrid=False)
    suptitle("12bp sample", fontsize=fontsize)
    gcf().subplots_adjust(top=0.9, bottom=0.1)
    gcf().subplots_adjust(hspace=0, wspace=0)
    for i,ax in enumerate(AX.ravel()):
        ax.set_xlabel(''); ax.set_ylabel(''); ax.set_title('')
        ax.set_yticks([])        
        ax.set_xticks([0,0.5,1])
        ax.text(0,3, "CH%d" % (i+1), fontsize=fontsize,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.text(0,1.5, "E = %.2f" % E_fit12[i], fontsize=fontsize)
                #bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.axvline(E_fit12[i], lw=3, color='k', ls='--', alpha=0.7)
    for ax in AX[-1,:]: ax.set_xlabel('E', fontsize=fontsize)
        
    
    
