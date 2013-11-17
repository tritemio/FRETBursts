# AX,s = dplot(dfs12, hist_fret, figsize=(5, 5.975),
# bins=r_[-0.2:1.2:0.05]+0.025, pgrid=False, name='')
ip.magic("run -i style.py")
suptitle("7bp sample", fontsize=fontsize)
gcf().subplots_adjust(top=0.9, bottom=0.1)
gcf().subplots_adjust(hspace=0., wspace=0.)
for i,ax in enumerate(AX.ravel()):
    ax.set_xlabel(''); ax.set_ylabel(''); ax.set_title('')
    #ax.set_yticks([])        
    #ax.set_xticks([0,0.5,1])
    ax.set_ylim(ymax=499)
    ax.text(0.8,0.8,"CH%d"%(i+1), fontsize=fontsize,transform = ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    #ax.text(0,1.5, "E = %.2f" % E_fit12[i], fontsize=fontsize)
    #        #bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    #ax.axvline(E_fit12[i], lw=3, color='k', ls='--', alpha=0.7)
    ax.figure.canvas.draw()
for ax in AX[-1,:]: ax.set_xlabel('Burst size', fontsize=fontsize)
    


