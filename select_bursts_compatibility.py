"""
This module is a compatibility layer for old scripts. It redefines all the
functions in select_bursts.py adding the prefix 'select_bursts_' to their name.

This module must be imported as:

    from select_bursts_compatibility import *

The new syntax to call a burst selection function is:

    select_bursts.nda()

This modules make the selection functions available with the old name:

    select_bursts_nda()

"""

import select_bursts
from utils.misc import pprint

def deprecate(function, old_name, new_name):
    def deprecated_function(*args, **kwargs):
        pprint("Function %s is deprecated, use %s instead.\n" %\
                (old_name, new_name))
        res = function(*args, **kwargs)
        return res
    return deprecated_function

for name, value in select_bursts.__dict__.items():
    if name == '__builtins__': continue
    if hasattr(select_bursts.__builtins__, name): continue
    if name.startswith('__'): continue
    print "Redefining `%s` as `%s`" % (name, 'select_bursts_'+name)
    globals()['select_bursts_'+name] = deprecate(value, name, 
                                                 'select_bursts_'+name)

