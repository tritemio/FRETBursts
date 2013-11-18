"""
This script run some consistency tests to check correcteness of processing.

To run the script you must load in a variable `d` a Data() object containing 
burst search results.
"""

# Test coherence between b_end() and b_iend()
for i in xrange(d.nch):
    assert (d.ph_times_m[i][b_iend(d.mburst[i])] == b_end(d.mburst[i])).all()

# Test for monotonic burst_start
for i in xrange(d.nch):
    assert (diff(b_start(d.mburst[i])) > 0).all()

# Test for monotonic burst_istart
for i in xrange(d.nch):
    assert (diff(b_istart(d.mburst[i])) > 0).all()

# Test for monotonic burst_end
for i in xrange(d.nch):
    assert (diff(b_end(d.mburst[i])) > 0).all()

## Test consistency between burst start, end and size
for mb in d.mburst:
    assert (mb[:,iiend] == mb[:,iistart]+mb[:,inum_ph]-1).all()
    assert (b_iend(mb) == b_istart(mb) + b_size(mb) - 1).all()


## Test that after fusing with ms=0 the sum of bursts sizes is that same
## as the number of ph in bursts (via burst selection)
if not hasattr(d, 'fuse'):
    #fused_mburst = mch_fuse_bursts(d.ph_times_m, d.mburst, ms=0)
    df = d.fuse_bursts(ms=0)
    for ph,mb in zip(d.ph_times_m, d.mburst):
        m = ph_select(ph, mb)
        assert m.sum() == b_size(mb).sum()










