#!/usr/bin/env python
"""
process_selection.py

Process 'Selection' tab for FlameTransfer app

Created December 2016 by C. Lapeyre (lapeyre@cerfacs.fr)
"""
import os

from os.path import isfile

from XDR2.exceptions import XDRException
from flametransferprocess import FlameTransferProcess
from process_library import process_library

def process_selection(pr):
    """Refresh the list of flames in the library"""
    ds = pr.ds

    try:
        selected_flame = ds.getSelection('selection', 'edit')[0]
    except IndexError:
        pr.log.warning("No action performed if no flame is checked!")
        ds.setValue('-', 'modified_flame')
        raise XDRException
    pr.log.debug('Selected flame: ' + selected_flame)
    if ds.getValue('action') == 'modify':
        ds.setValue(selected_flame, 'modified_flame')
        if selected_flame == 'new':
            pr.log.info("Creating new flame")
        else:
            metas = pr.get_metas()[selected_flame]
            sh_args = metas["shape_params"].split()
            ds.setValue(to_ds_list('x y z'.split(), metas["ref_point"].split()), "ptref_list")
            ds.setValue(to_ds_list('x y z'.split(), metas["ref_vect"].split()), "vecref_list")
            ds.removeNode("xor_ndim")
            ndim = "three_d" if "3D" in metas["generation_method"] else "two_d"
            ds.addChild("xor_ndim", ndim, "modify")
            ds.addChild("xor_flame_geo", metas["generation_method"], "xor_ndim")
            ds.addChild(metas["generation_method"], "", "xor_flame_geo")
            if metas["generation_method"] == "analytical2D_parallelogram":
                ds.addChild("ana_flame_center", to_ds_list('x y'.split(), sh_args[0:2]), metas["generation_method"])
                ds.addChild("ana_flame_vector1", to_ds_list('u1_x u1_y'.split(), sh_args[2:4]), metas["generation_method"])
                ds.addChild("ana_flame_vector2", to_ds_list('u2_x u2_y'.split(), sh_args[4:]), metas["generation_method"])
            elif metas["generation_method"] == "analytical2D_circle":
                ds.addChild("ana_flame_center", to_ds_list('x y'.split(), sh_args[0:2]), metas["generation_method"])
                ds.addChild("ana_flame_radius", str(sh_args[2]), metas["generation_method"])
            elif metas["generation_method"] == "analytical3D_sphere":
                ds.addChild("ana_flame_center", to_ds_list('x y z'.split(), sh_args[0:3]), metas["generation_method"])
                ds.addChild("ana_flame_radius", str(sh_args[3]), metas["generation_method"])
            elif metas["generation_method"] == "analytical3D_cylinder":
                ds.addChild("ana_flame_center", to_ds_list('x y z'.split(), sh_args[0:3]), metas["generation_method"])
                ds.addChild("ana_flame_radius", str(sh_args[3]), metas["generation_method"])
                ds.addChild("ana_flame_vector", to_ds_list('u_x u_y u_z'.split(), sh_args[4:]), metas["generation_method"])
            else:
                pr.log.error("Unknown generation method {}. Cannot set default values in 'modify' tab".format(metas["generation_method"]))
                raise ValueError
    elif selected_flame == "new":
        pr.log.warning("'new' flame only supports 'Modify'. It is then filled in the 'modify' tab")
    else:
        if ds.getValue('action') == 'duplicate':
            pr.duplicate_flame(selected_flame, selected_flame + '_duplicate')
        elif ds.getValue('action') == 'delete':
            pr.def_flame(selected_flame)


def to_ds_list(keys, values):
    out = []
    for k, v in zip(keys, values):
        out.extend([str(k), str(v)])
    return ";".join(out)


if __name__ == '__main__':
    pr = FlameTransferProcess('dataset.xml')
    process_selection(pr)
    process_library(pr)
    pr.finish()
