"""
This version of burst search maps exactly to Labview-FPGA constructs.

Used for logical validation of the FPGA algorithm.
"""


m = 3

prefilling = True   # (bool) Non-default init
buffered_ph = 0     # (uint8)
in_burst = False    # (bool)        
burst_start = 0     # (uint32)
num_ph = 0          # (uint16)
num_ph_donor = 0    # (uint16)
donor_fifo = []
t_start_fifo = []

def donor_counter(donor_fifo):
    return reduce(lambda x,y:x+y, donor_fifo)

def fpga_burst_search(data_avail, t, donor, m, T, L):
    """Function executed at every clock cycle."""
    global prefilling, buffered_ph, in_burst, burst_start, num_ph, num_ph_donor
    global t_start, t_secondlast

    burst_available = False
    out_burst_start,out_burst_width,out_burst_size,out_burst_size_donor=[0,0,0,0]
    if data_avail: # T
        t_start_fifo.append(t)
        donor_fifo.append(donor)
        if prefilling: # T T
            buffered_ph += 1
            prefilling = (buffered_ph < m-1)
        else: # T F
            t_start = t_start_fifo.pop(0)
            donor_fifo.pop(0)
            if t-t_start <= T: # T F T      (ABOVE the min rate)
                if not in_burst: # T F T F
                    in_burst = True
                    burst_start = t_start
                    num_ph = m
                    num_ph_donor = donor_counter(donor_fifo)
                else: # T F T T
                    num_ph += 1
                    if donor: num_ph_donor += 1
            else: # T F F                   (BELOW the min rate)
                if in_burst: # T F F T
                    in_burst = False
                    if num_ph >= L: # T F F T T
                        burst_available = True
                        out_burst_start = burst_start
                        burst_end = t_secondlast
                        out_burst_width = burst_end - burst_start
                        out_burst_size = num_ph
                        out_burst_size_donor = num_ph_donor
                    else: # T F F T F
                        pass
                else: # T F F F
                    pass
                num_ph = 0
                num_ph_donor = 0
     
    t_secondlast = t
    return burst_available, \
        out_burst_start, out_burst_width, out_burst_size, out_burst_size_donor

def test():
    _ip.magic("run -i burst")
    ph0 = ph_times_m[0]
    don0 = -red[0]
    burst_fpga = []
    pprint('\n - FPGA burst search ... ')
    for i in range(ph0.size):
        l = fpga_burst_search(True, ph0[i], don0[i], 3, 10000, 20)
        if l[0]: burst_fpga.append(l[1:])
    pprint('[DONE]\n\n')

    m0 = array(burst_fpga)
    m1 = ba(ph0, 20, 3, 10000)
    
    assert m1.shape == m0.shape
    assert (m1[:,0] == m0[:,0]).sum() == m1.shape[0]
    assert (m1[:,1] == m0[:,1]).sum() == m1.shape[0]
    assert (m1[:,2] == m0[:,2]).sum() == m1.shape[0]
    print "\n Test passed: FPGA burst search gives the same results!"

