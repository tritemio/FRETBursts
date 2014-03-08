#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
SPC Format (Beker & Hickl)
--------------------------

48-bit element in little endian (<) format

Drawing (note: each char represents 2 bits)::

    bit: 64        48                          0
         0000 0000 XXXX XXXX XXXX XXXX XXXX XXXX
                   '-------' '-------' '-------'
    uint16:         data[2]   data[1]   data[0]

         0000 0000 XXXX XXXX XXXX XXXX XXXX XXXX
                   '-------' '--' '--'   '-----'
                       a      c    b        d

    macrotime = [ b  ] [     a     ]  (24 bit)
    detector  = [ c  ]                (8 bit)
    nanotime  = [  d  ]               (12 bit)

    overflow bit: 13, bit_mask = 2^(13-1) = 4096
"""

from numpy import *

def load_spc(fname, return_extra_data=False):
    f = open(fname, 'rb')
    raw_data = f.read()

    ## Each element is 48bit, which will be imported as 3 uint16
    bytes_per_element = 6
    N_elements = len(raw_data)/bytes_per_element
    data = ndarray(shape=(3*N_elements,), buffer=raw_data,dtype='<u2')
    data = data.reshape(data.size/3,3) # every row is a 48bit element

    b = bitwise_and(data[:,1], int('0xFF', 16))
    ph_times = left_shift(b.astype(uint64),16) + data[:,2] # [ b ][   a   ]
    
    det = right_shift(data[:,1],8)
    nanotime = 4095 - bitwise_and(data[:,0], int('0xFFF',16))
    extra_data = right_shift(data[:,0],12).astype(uint8)
    
    # extract the 13-th bit of every elemente in data[:,0]
    overflow = bitwise_and(right_shift(data[:,0],13), 1)
    overflow = cumsum(overflow, dtype=uint64)

    # concatenate of the MSB given by the overflow integration 
    ph_times += left_shift(overflow,24)
    
    # sanity checks
    assert ph_times.shape == det.shape
    assert ph_times.shape == nanotime.shape

    # The first ph has always detector==1, seems odd so I exclude it
    if return_extra_data:
        return ph_times[1:], det[1:], nanotime[1:], extra_data[1:]
    else:    
        return ph_times[1:], det[1:], nanotime[1:]


