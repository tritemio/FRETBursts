"""
Utility functions
"""

import sys

def clk_to_s(t_ck, clk_p=12.5*1e-9):
    """Convert clock cycles to seconds."""
    return t_ck*clk_p

def pprint(s):
    """Print immediately, even if inside a busy loop."""
    sys.stdout.write(s)
    sys.stdout.flush()


