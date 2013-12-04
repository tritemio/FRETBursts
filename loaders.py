"""
Functions to load supported file formats into a Data() object.

These are high-level helper functions that just pack the data in a Data()
object. The low-level format decoding functions are in dataload folder.
"""

import cPickle as pickle
from dataload.multi_ch_reader import load_data_ordered16
from dataload.smreader import load_sm


##
# Multi-spot loader functions
#
def load_multispot8(fname, bytes_to_read=-1, swap_D_A=True, BT=0, gamma=1.):
    """Load a 8-ch multispot file and return a Data() object. Cached version.
    """
    fname_c = fname + '_cache.pickle'
    try:
        var = pickle.load(open(fname_c, 'rb'))
        dx = Data(fname=fname, clk_p=12.5e-9, nch=8, BT=BT, gamma=gamma)
        dx.add(ph_times_m=var['ph_times_m'], A_em=var['A_em'], ALEX=False)
        pprint(" - File loaded from cache: %s\n" % fname)
    except IOError:
        dx = load_multispot8_core(fname, bytes_to_read=bytes_to_read,
                                  swap_D_A=swap_D_A, BT=BT, gamma=gamma)
        D = {'ph_times_m': dx.ph_times_m, 'A_em': dx.A_em}
        pprint(" - Pickling data ... ")
        pickle.dump(D, open(fname_c, 'wb'), -1)
        pprint("DONE\n")
    return dx

def load_multispot8_core(fname, bytes_to_read=-1, swap_D_A=True, BT=0,
                         gamma=1.):
    """Load a 8-ch multispot file and return a Data() object.
    """
    dx = Data(fname=fname, clk_p=12.5e-9, nch=8, BT=BT, gamma=gamma)
    ph_times_m, A_em, ph_times_det = load_data_ordered16(fname=fname,
            n_bytes_to_read=bytes_to_read, swap_D_A=swap_D_A)
    dx.add(ph_times_m=ph_times_m, A_em=A_em, ALEX=False)
    return dx

##
# usALEX loader functions
#

# Build masks for the alternating periods
def select_outer_range(times, period, edges):
    return ((times%period) > edges[0]) + ((times%period) < edges[1])

def select_inner_range(times, period, edges):
    return ((times%period) > edges[0]) * ((times%period) < edges[1])

def load_usalex(fname, bytes_to_read=-1, swap_D_A=True, BT=0,
                         gamma=1., header=166):
    """Load a usALEX file and return a Data() object.
    """
             
    print " - Loading '%s' ... " % fname
    ph_times_t, det_t = load_sm(fname, header=header) 
    print " [DONE]\n"
    
    donor_ch, accept_ch = 1, 0
    if swap_D_A: 
        print '- Swaping D and A'
        D_ch, A_ch = accept_ch, donor_ch
    else:
        D_ch, A_ch = donor_ch, accept_ch
    
    DONOR_ON = (2850, 580)
    ACCEPT_ON = (930, 2580)
    switch_window = 2000

    dx = Data(fname=fname, clk_p=12.5e-9, nch=1, BT=BT, gamma=gamma, 
              ALEX=True,
              D_ON=DONOR_ON, A_ON=ACCEPT_ON, switch_window=switch_window,
              ph_times_t=[ph_times_t], D_em_t=[(det_t == D_ch)],
              )
    return dx

def usalex_apply_period(d):  
    
    ph_times_t = d.ph_times_t[0]
    period = 2*d.switch_window
    
    # Create masks for donor and acceptor channels
    donor_t_mask = d.D_em_t[0]
    accept_t_mask = -d.D_em_t[0]
    print "#donor: %d  #acceptor: %d \n" % (donor_t_mask.sum(), 
            accept_t_mask.sum())

    # Build masks for excitation windows
    if d.D_ON[0] < d.D_ON[1]: 
        donor_ex_t = select_inner_range(ph_times_t, period, d.D_ON)
    else:
        donor_ex_t = select_outer_range(ph_times_t, period, d.D_ON)
    if d.A_ON[0] < d.A_ON[1]:
        accept_ex_t = select_inner_range(ph_times_t, period, d.A_ON)
    else:
        accept_ex_t = select_outer_range(ph_times_t, period, d.A_ON)

    # Safety check: each ph is either D or A ex (not both)
    assert (donor_ex_t*accept_ex_t == False).any() 
    
    # (3) Burst search on donor_excitation+acceptor_excitation
    mask = donor_ex_t + accept_ex_t
    
    #ph_burst_search = "donor_accept_ex"

    # Assign the new ph selection mask
    ph_times = ph_times_t[mask]
    d_em = donor_t_mask[mask]
    a_em = accept_t_mask[mask]
    d_ex = donor_ex_t[mask]
    a_ex = accept_ex_t[mask]
    
    assert d_em.sum() + a_em.sum() == ph_times.size
    assert (d_em * a_em).any() == False
    assert a_ex.size == a_em.size == d_ex.size == d_em.size == ph_times.size
    
    ph_times_m, A_ex, D_ex, A_em, D_em = [ph_times],[a_ex],[d_ex],[a_em],[d_em]
    d.add(ph_times_m=ph_times_m, 
          D_em=D_em, A_em=A_em, D_ex=D_ex, A_ex=A_ex,
          )
    assert d.ph_times_m[0].size == d.A_em[0].size
    assert d.switch_window*2 == period
    
    return d
              


