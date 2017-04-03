#!/usr/bin/env python
"""
process_library.py

Process 'Library' tab (actually 'Refresh list' button) for FlameTransfer app

Created December 2016 by C. Lapeyre (lapeyre@cerfacs.fr)
"""
from flametransferprocess import FlameTransferProcess


def process_library(pro):
    """Refresh the list of flames in the library"""
    dst = pro.ds
    action = dst.getValue("xor_lib_action")
    if action == "refresh_flame_list":
        pass  # Refresh done on process initialization
    if action == "create_new_flame":
        name = dst.getValue("input_flame_name", "xor_lib_action")
        pro.add_libobj(name)
        dst.setValue(name, "cho_write_to_flame")
    if action == "import_flame":
        path = dst.getValue("inputfile_import_flame", "xor_lib_action")
        pro.import_libobj(path)
        dst.setValue(pro.libobjs[-1], "cho_write_to_flame")
    if action == "replicate_flame":
        name = dst.getValue("cho_flame", "xor_lib_action")
        typ = dst.getValue("xor_transform")
        if typ == "none":
            new_name = (
                    dst.getValue("input_flame_name") if
                    dst.getValue("flag_rename_duplicate")
                    else None)
            pro.duplicate_libobj(name, new_name)
        elif typ == "translate":
            vector = [float(e) for e in
                      dst.getValue("list_translate_flame").split(';')[1::2]]
            final_nb = int(dst.getValue("int_nb_flames_final"))
            pro.replicate_translate(name, vector, final_nb)
        elif typ == "rotate":
            angle = float(dst.getValue("dbl_rotate_angle"))
            final_nb = int(dst.getValue("int_nb_flames_final"))
            pro.replicate_rotate(name, angle, final_nb)
    if action == "rename_flame":
        src = dst.getValue("cho_flame", "xor_lib_action")
        dest = dst.getValue("input_flame_name", "xor_lib_action")
        pro.rename_libobj(src, dest)
        dst.setValue(dest, "cho_write_to_flame")
    if action == "delete_flame":
        name = dst.getValue("cho_flame", "xor_lib_action")
        pro.del_libobj(name)
    pro.update_libobjs()


if __name__ == '__main__':
    pro = FlameTransferProcess('dataset.xml')
    process_library(pro)
    pro.finish()
