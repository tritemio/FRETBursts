"""
Functions to load supported file formats into a Data() object.

These are high-level helper functions that just pack the data in a Data()
object. The low-level format decoding functions are in dataload folder.
"""

import cPickle as pickle
from dataload.multi_ch_reader import load_data_ordered16


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


