#!/usr/bin/env python
"""
process_library.py

Process 'Library' tab (actually 'Refresh list' button) for FlameTransfer app

Created December 2016 by C. Lapeyre (lapeyre@cerfacs.fr)
"""
from flametransferprocess import FlameTransferProcess


def process_library(pr):
    """Refresh the list of flames in the library"""
    ds = pr.ds
    action = ds.getValue("xor_lib_action")
    if action == "refresh_flame_list":
        pass  # Refresh done on process initialization
    if action == "create_new_flame":
        name = ds.getValue("input_flame_name", "xor_lib_action")
        pr.add_libobj(name)
        ds.setValue(name, "cho_write_to_flame")
    if action == "import_flame":
        path = ds.getValue("inputfile_import_flame", "xor_lib_action")
        pr.import_libobj(path)
        ds.setValue(pr.libobjs[-1], "cho_write_to_flame")
    if action == "replicate_flame":
        name = ds.getValue("cho_flame", "xor_lib_action")
        typ = ds.getValue("xor_transform")
        if typ == "none":
            new_name = (
                    ds.getValue("input_flame_name") if
                    ds.getValue("flag_rename_duplicate")
                    else None)
            pr.duplicate_libobj(name, new_name)
        elif typ == "translate":
            vector = [float(e) for e in
                      ds.getValue("list_translate_flame").split(';')[1::2]]
            final_nb = int(ds.getValue("int_nb_flames_final"))
            pr.replicate_translate(name, vector, final_nb)
        elif typ == "rotate":
            angle = float(ds.getValue("dbl_rotate_angle"))
            final_nb = int(ds.getValue("int_nb_flames_final"))
            pr.replicate_rotate(name, angle, final_nb)
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
