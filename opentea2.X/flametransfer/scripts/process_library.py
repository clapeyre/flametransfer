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
        pr.add_libobj(name)
        ds.setValue(name, "cho_write_to_flame")
    if action == "import_flame":
        path = ds.getValue("inputfile_import_flame", "xor_lib_action")
        pr.import_libobj(path)
        ds.setValue(pr.libobjs[-1], "cho_write_to_flame")
    if action == "duplicate_flame":
        name = ds.getValue("cho_flame", "xor_lib_action")
        pr.duplicate_libobj(name)
        ds.setValue(name+'_duplicate', "cho_write_to_flame")
    if action == "rename_flame":
        src = ds.getValue("cho_flame", "xor_lib_action")
        dest = ds.getValue("input_flame_name", "xor_lib_action")
        pr.rename_libobj(src, dest)
        ds.setValue(dest, "cho_write_to_flame")
    if action == "delete_flame":
        name = ds.getValue("cho_flame", "xor_lib_action")
        pr.del_libobj(name)
    pr.update_libobjs()

if __name__ == '__main__':
    pr = FlameTransferProcess('dataset.xml')
    process_library(pr)
    pr.finish()
