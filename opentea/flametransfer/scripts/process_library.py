#!/usr/bin/env python
"""
process_library.py

Process 'Library' tab (actually 'Refresh list' button) for FlameTransfer app

Created December 2016 by C. Lapeyre (lapeyre@cerfacs.fr)
"""
from glob import glob

from XDR2 import XDRUnknownValue
from flametransferprocess import FlameTransferProcess

def process_library(pr):
    """Refresh the list of flames in the library"""
    ds = pr.ds
    action = ds.getValue("xor_lib_action")
    if action == "refresh_flame_list":
        pass # This is done on process initialization
    if action == "create_new_flame":
        name = ds.getValue("input_flame_name", "xor_lib_action")
        pr.add_flame(name)
    if action == "import_flame":
        path = ds.getValue("inputfile_import_flame", "xor_lib_action")
        pr.import_flame(path)
    if action == "duplicate_flame":
        name = ds.getValue("cho_flame", "xor_lib_action")
        pr.duplicate_flame(name)
    if action == "rename_flame":
        src = ds.getValue("cho_flame", "xor_lib_action")
        dest = ds.getValue("input_flame_name", "xor_lib_action")
        pr.rename_flame(src, dest)
    if action == "delete_flame":
        name = ds.getValue("cho_flame", "xor_lib_action")
        pr.del_flame(name)
    pr.update_flames()

if __name__ == '__main__':
    pr = FlameTransferProcess('dataset.xml')
    process_library(pr)
    pr.finish()
