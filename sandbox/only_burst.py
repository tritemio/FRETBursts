"""
Script used in older version for part of the analysis. Now obsolete.
"""

#
# Script called after burst.py to do the burst search on the loaded data
#

## Test the burst fusion
#mburstm = [d.fuse_bursts(i) for i in range(d.nch)]
#bmask = hstack([(burst_separation(d,ich=1) <= 0), (False,)]

## Burst search parameters
#T_us = 1000; m = 3;
#L = 20

# Save everything in a variable of type Data (see burstlib.py)
#d = Data(fname=short_fname, nch=4, clk_p=clk_p, BT=BT,  
#        ALEX=False, ph_times=ph_times, det=det, 
#        ph_times_m=ph_times_m,red=red,
#        ratem=ratem, rates=rates, 
#        )

## Burst search on 4 merged channels
#mburst, T_ch = seek_burst_ap(ph_times_m, rates4, m, L, P_user=0.1)
#Num_green, Num_red = mch_count_burst_donor_accept(mburst, red)

# Put all the data in an object
#d = Data(fname=short_fname, nch=4, clk_p=clk_p, BT=BT, corrected=False,
#        #ph_times=ph_times, det=detector, masks=masks,
#        ALEX=False, ph_times_m=ph_times_m, red=red, 
#        ratem=rates4, rates=rates8, 
#        L=L, m=m, T=T_ch, 
#        mburst=mburst, nd=Num_green, na=Num_red)

# Compute corrections
#d.cal_corrections()
#background_subtract(d)
#bleed_through_correct(d, BT)

# Compute FRET efficiency
#d.cal_fret()

# Save some statistics
name = short_fname.replace('/','_')[:-4]
s = print_burst_stats(d)
try:
    # May fail if there is no minimum between donor and acceptor population
    # in the FRET histogram
    s += print_fret_stats(d) 
except:
    pass
#log_dir = '../log/' -> Now defined in path_def.py
print >> open(log_dir+name+'_L'+str(L)+'.txt','w'), s
print s

#fig_dir = '../figure/' -> Now defined in path_def.py
if not os.path.exists(fig_dir+name):
    pprint('\n - Creating the figure folder "%s" ...' % (fig_dir+name))
    os.mkdir(fig_dir+name)
    pprint('[DONE]\n')
else:
    pprint('\n - Figure folder already present.\n')
#multi_plot(d, save=False, subdir=name)

#savefig(fig_dir+name+'/'+name+'_'+ptype+'_L'+str(L)+'.png')

