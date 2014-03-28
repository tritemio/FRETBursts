#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
SM Format, written by the LV program in WeissLab us-ALEX setup
--------------------------------------------------------------

A SM file is composed by two parts:

  - 166 bytes: an header
  - the remainig bytes: the data
  - 26 bytes of trailing cruft

The data is a list of 96-bit records that contain for each ph (count)
the detection time and the detector.

The first 64-bit of each record are the ph time, and the remainig 16-bit
are the detector.

::
              64-bit ph time        detector
            -------------------     -------
           |                   |   |       |
           XXXX XXXX - XXXX XXXX - XXXX XXXX
     bit:  0      32   0      32   0      32
           '-------'   '-------'   '-------'
            data[0]     data[1]     data[2]

The data is in unsigned big endian (>) format.

For efficiency the data is read as a list of B.E. (>) uint32 and then
recomposed in ph_times and det. (NB: Old not very efficient approach)

The proper way to read the data is to read the byte-stream and interpret
it as a record array in which each element is 12 bytes.
"""

import numpy as np


def load_sm(fname, header=166):
    try: f = open(fname, 'rb')
    except IOError: f = open(fname+'.sm', 'rb')
    head = f.read(header)
    raw_data = f.read()
    raw_data = raw_data[:4*(len(raw_data)/4)]
    data = np.fromstring(raw_data,'>u4')
    data = data.reshape(data.size/3,3)[:-20,:]
    #ph_times = left_shift(data[:,0].astype(uint64),32)+data[:,1].astype(uint64)
    ph_times = np.left_shift(data[:,0].astype('int64'),32) + \
                data[:,1].astype('int64')
    det = data[:,2].astype('uint8') # 8 bit are enough (also bool!)
    return ph_times, det

def load_sm_new(fname, header=166):
    try: f = open(fname, 'rb')
    except IOError: f = open(fname+'.sm', 'rb')
    f.seek(header)
    byte_string = f.read()

    # Description of the record element in the file
    sm_dtype = np.dtype([('timestamp', '>i8'), ('detector', '>u4')])

    # Remove the end of the file
    end_field1 = 4
    end_str = 'End Of Run'
    end_field2 = 12
    valid_size = len(byte_string) - end_field1 - len(end_str) - end_field2

    # View of the binary dtaa as an array (no copy performed)
    data = np.frombuffer(byte_string[:valid_size], dtype=sm_dtype)
    return data['timestamp'], data['detector']

def _decode_header(header):
    """Decode the header of a .sm file. UNUSED."""
    # List of ASCII strings in the header
    str_list = [
            'Simple', 'Arrival Time Counter',
            'Time High', 'Time Low',
            'Ch 1', 'Ch 2'
            ]

    # List of ASCII byte fields lengths
    field_list = [12, 8, 12, 24, 55, 4]

    str_len = len(''.join(str_list))
    field_len = reduce(lambda x,y: x+y, field_list)

    assert len(str_list) == len(field_list)
    assert str_len + field_len == len(header)

    Fields = []
    Strings = []
    pos = 0
    for field, string in zip(field_list, str_list):
        # Read the sequence of byte fields and strings
        d = np.frombuffer(header[pos:pos+field], count=field/4, dtype='>u4')
        pos += field
        s = header[pos:pos+len(string)]
        pos += len(string)
        Fields.append(d)
        Strings.append(s)

    return Strings, Fields


