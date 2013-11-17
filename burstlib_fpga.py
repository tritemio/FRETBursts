##
# FPGA Functions
#
def write_fpga_file(d, fname):
    counter = 0
    out_file = open(fname, 'wb')
    header = 'channels = 8\nwords per burst = 6\nclock frequency (MHz) = 80\n'
    out_file.write(header)
    for ch in range(d.nch):
        N = d.mburst[ch].shape[0]
        col1 = zeros(N, dtype=uint32)
        for i in xrange(N):
            col1[i] = bitwise_or(uint32(ch), left_shift(uint8(counter),16))
            col1[i] = bitwise_or(left_shift(0xF0, 24), col1[i])
            counter += 1
        bursts = d.mburst[ch][:,:3].astype(uint32)
        col6 = uint32(d.T[0]/d.clk_p) * ones(N, dtype=uint32)
        data = vstack([col1, bursts.T, d.nd[ch].astype(uint32),
            col6]).T
        s = data.tostring()
        print len(s)
        out_file.write(s)
    out_file.close()
    return data, s

#def load_fpga_file(fname):
#    f = open(fname)
#    f.readline()
#    f.readline()
#    f.readline()
#    s = f.read()
#    N = len(s)/4/6
#    dat = fromstring(s, dtype='<u4')
#    return dat.reshape(dat.size/6,6)
