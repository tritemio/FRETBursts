"""
Functions to load supported file formats into a Data() object.

These are high-level helper functions that just pack the data in a Data()
object. The low-level format decoding functions are in dataload folder.
"""

import os
import numpy as np
import cPickle as pickle
from dataload.multi_ch_reader import load_data_ordered16
from dataload.smreader import load_sm
from dataload.manta_reader import (load_manta_timestamps,
                                   load_xavier_manta_data,
                                   get_timestamps_detectors,
                                   #process_timestamps,
                                   process_store,)
from dataload.pytables_array_list import PyTablesList

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

def load_multispot48_simple(fname, BT=0, gamma=1.,
                     i_start=0, i_stop=None, debug=False):
    """Load a 48-ch multispot file and return a Data() object.
    """
    dx = Data(fname=fname, clk_p=10e-9, nch=48, BT=BT, gamma=gamma)
    ph_times_m, big_fifo, ch_fifo = load_manta_timestamps(
                fname, i_start=i_start, i_stop=i_stop, debug=debug)
    A_em = [True] * len(ph_times_m)
    dx.add(ph_times_m=ph_times_m, A_em=A_em, ALEX=False)
    big_fifo_full = np.array([b.any() for b in big_fifo]).any()
    ch_fifo_full = np.array([b.any() for b in ch_fifo]).any()
    if big_fifo_full:
        print 'WARNING: Big-FIFO full, flags saved in Data()'
        dx.add(big_fifo=big_fifo)
    if ch_fifo_full:
        print 'WARNING: CH-FIFO full, flags saved in Data()'
        dx.add(ch_fifo=ch_fifo)
    return dx

def load_multispot48(fname, BT=0, gamma=1.,
                     i_start=0, i_stop=None, debug=False):
    """Load a 48-ch multispot file and return a Data() object.
    """
    fname_h5 = fname if fname.endswith('hdf5') else fname[:-3]+'hdf5'

    if not os.path.exists(fname):
        raise IOError('Data file "%s" not found' % fname)

    if os.path.exists(fname_h5):
        ## There is a HDF5 file
        pprint(' - Loading HDF5 file: %s ... ' % fname_h5)
        ph_times_m = PyTablesList(fname_h5, group_name='timestamps_list')

        big_fifo = PyTablesList(ph_times_m.data_file,
                                group_name='big_fifo_full_list',
                                parent_node='/timestamps_list')

        ch_fifo = PyTablesList(ph_times_m.data_file,
                                group_name='small_fifo_full_list',
                                parent_node='/timestamps_list')
        pprint('DONE.\n')
    else:
        pprint(' - Loading file: %s ... ' % fname)
        ## Load data from raw file and store it in a HDF5 file
        data = load_xavier_manta_data(fname, i_start=i_start, i_stop=i_stop, 
                                      debug=debug)
        pprint('DONE.\n - Extracting timestamps and detectors ... ')
        timestamps, det = get_timestamps_detectors(data, nbits=24)
        pprint('DONE.\n - Processing and storing ... ')
        ph_times_m, big_fifo, ch_fifo = process_store(timestamps, det,
                        out_fname=fname_h5, fifo_flag=True, debug=False)
        pprint('DONE.\n')
    ## Current data has only acceptor ch
    A_em = [True] * len(ph_times_m)

    dx = Data(fname=fname, clk_p=10e-9, nch=48, BT=BT, gamma=gamma)
    dx.add(ph_times_m=ph_times_m, A_em=A_em, ALEX=False)
    big_fifo_full = np.array([b[:].any() for b in big_fifo]).any()
    ch_fifo_full = np.array([b[:].any() for b in ch_fifo]).any()
    if big_fifo_full:
        print 'WARNING: Big-FIFO full, flags saved in Data()'
        dx.add(big_fifo=big_fifo)
    if ch_fifo_full:
        print 'WARNING: CH-FIFO full, flags saved in Data()'
        dx.add(ch_fifo=ch_fifo)
    return dx


##
# usALEX loader functions
#

# Build masks for the alternating periods
def _select_outer_range(times, period, edges):
    return ((times % period) > edges[0]) + ((times % period) < edges[1])

def _select_inner_range(times, period, edges):
    return ((times % period) > edges[0]) * ((times % period) < edges[1])

def _select_range(times, period, edges):
    return _select_inner_range(times, period, edges) if edges[0] < edges[1] \
            else _select_outer_range(times, period, edges)

def load_usalex(fname, BT=0, gamma=1., header=166, bytes_to_read=-1):
    """Load a usALEX file and return a Data() object.

    To load usALEX data follow this pattern:

        d = load_usalex(fname=fname, BT=0, gamma=1.)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580), alex_period=4000)
        plot_alternation_hist(d)

    If the plot looks good apply the alternation with:

        usalex_apply_period(d)
    """
    print " - Loading '%s' ... " % fname
    ph_times_t, det_t = load_sm(fname, header=header)
    print " [DONE]\n"

    DONOR_ON = (2850, 580)
    ACCEPT_ON = (930, 2580)
    alex_period = 4000

    dx = Data(fname=fname, clk_p=12.5e-9, nch=1, BT=BT, gamma=gamma,
              ALEX=True,
              D_ON=DONOR_ON, A_ON=ACCEPT_ON, alex_period=alex_period,
              ph_times_t=ph_times_t, det_t=det_t, det_donor_accept=(0, 1),
              )
    return dx

def usalex_apply_period(d, delete_ph_t=True, remove_d_em_a_ex=False):
    """Applies the alternation period previously set.

    To load usALEX data follow this pattern:

        d = load_usalex(fname=fname, BT=0, gamma=1.)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580), alex_period=4000)
        plot_alternation_hist(d)

    If the plot looks good apply the alternation with:

        usalex_apply_period(d)
    """
    donor_ch, accept_ch  = d.det_donor_accept
    # Remove eventual ch different from donor or acceptor
    d_ch_mask_t = (d.det_t == donor_ch)
    a_ch_mask_t = (d.det_t == accept_ch)
    valid_mask = d_ch_mask_t + a_ch_mask_t
    ph_times_val = d.ph_times_t[valid_mask]
    d_ch_mask_val = d_ch_mask_t[valid_mask]
    a_ch_mask_val = a_ch_mask_t[valid_mask]
    assert (d_ch_mask_val + a_ch_mask_val).all()
    assert not (d_ch_mask_val * a_ch_mask_val).any()

    # Build masks for excitation windows
    d_ex_mask_val = _select_range(ph_times_val, d.alex_period, d.D_ON)
    a_ex_mask_val = _select_range(ph_times_val, d.alex_period, d.A_ON)
    # Safety check: each ph is either D or A ex (not both)
    assert not (d_ex_mask_val * a_ex_mask_val).any()

    mask = d_ex_mask_val + a_ex_mask_val  # Removes alternation transients

    # Assign the new ph selection mask
    ph_times = ph_times_val[mask]
    d_em = d_ch_mask_val[mask]
    a_em = a_ch_mask_val[mask]
    d_ex = d_ex_mask_val[mask]
    a_ex = a_ex_mask_val[mask]

    if remove_d_em_a_ex:
        # Removes donor-ch photons during acceptor excitation
        mask = a_em + d_em*d_ex
        assert (mask == -(a_ex*d_em)).all()

        ph_times = ph_times[mask]
        d_em = d_em[mask]
        a_em = a_em[mask]
        d_ex = d_ex[mask]
        a_ex = a_ex[mask]

    assert d_em.sum() + a_em.sum() == ph_times.size
    assert (d_em * a_em).any() == False
    assert a_ex.size == a_em.size == d_ex.size == d_em.size == ph_times.size
    print "#donor: %d  #acceptor: %d \n" % (d_em.sum(), a_em.sum())

    d.add(ph_times_m=[ph_times],
          D_em=[d_em], A_em=[a_em], D_ex=[d_ex], A_ex=[a_ex],)

    assert d.ph_times_m[0].size == d.A_em[0].size

    if delete_ph_t:
        d.delete('ph_times_t')
        d.delete('det_t')
    return d



