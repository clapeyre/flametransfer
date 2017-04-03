#!/usr/bin/env python
"""
process_modify.py

Process 'Modify' tab for FlameTransfer app

Created December 2016 by C. Lapeyre (lapeyre@cerfacs.fr)
"""

import json
import numpy as np

from os.path import isfile, basename

from XDR2 import XDRInterrupt
from flametransferprocess import FlameTransferProcess

def process_modify(pr):
    """Modify existing or new virgin flame"""
    ds = pr.ds

    label = ds.getValue("cho_write_to_flame")
    put_files = []
    flame_geo = ds.getValue("xor_flame_geo", "modify")
    if flame_geo == "analytical2D_disc":
        center = ds.getListValue("ana_flame_center", "modify")[1::2]
        radius = ds.getValue("ana_flame_radius", "modify")
        script = ["ge disc " + label,
                  "{0[0]} {0[1]}".format(center),
                  "{}".format(radius)]
    elif flame_geo == "analytical2D_rectangle":
        pt_min = ds.getListValue("ana_flame_pt_min", "modify")[1::2]
        pt_max = ds.getListValue("ana_flame_pt_max", "modify")[1::2]
        script = ["ge rect " + label]
        script += ["{0[0]} {0[1]}".format(vect)
                   for vect in [pt_min, pt_max]]
    elif flame_geo == "analytical3D_brick":
        pt_min = ds.getListValue("ana_flame_pt_min", "modify")[1::2]
        pt_max = ds.getListValue("ana_flame_pt_max", "modify")[1::2]
        script = ["ge brick " + label]
        script += ["{0[0]} {0[1]} {0[2]}".format(vect)
                   for vect in [pt_min, pt_max]]
    elif flame_geo == "analytical3D_sphere":
        center = ds.getListValue("ana_flame_center", "modify")[1::2]
        radius = ds.getValue("ana_flame_radius", "modify")
        script = ["ge sphere " + label,
                  "{0[0]} {0[1]} {0[2]}".format(center),
                  "{}".format(radius)]
    elif flame_geo == "analytical3D_cylinder":
        center = ds.getListValue("ana_flame_center", "modify")[1::2]
        radius = ds.getValue("ana_flame_radius", "modify")
        vector = ds.getListValue("ana_flame_vector", "modify")[1::2]
        pr.log.debug("Creating cylinder with center " + repr(center))
        pr.log.debug("                       radius " + repr(radius))
        pr.log.debug("                       vector " + repr(vector))
        script = ["ge cylinder " + label,
                  "{0[0]} {0[1]} {0[2]}".format(center),
                  "{}".format(radius),
                  "{0[0]} {0[1]} {0[2]}".format(vector)]
    elif flame_geo in ["avbp_scalar_threshold_3D", "avbp_scalar_threshold_2D"]:
        avbp_sol = ds.getValue("avbp_sol", "xor_flame_geo", "modify")
        avbp_mesh = ds.getValue("avbp_mesh", "xor_flame_geo", "modify")
        scal = ds.getValue("scal", "xor_flame_geo", "modify")
        thresh = ds.getValue("thresh", "xor_flame_geo", "modify")
        script = ["ge avbp " + label,
                  basename(avbp_mesh),
                  basename(avbp_sol),
                  scal,
                  thresh,
                  "set static avbp_mesh string",
                  avbp_mesh,
                  "set static avbp_sol string",
                  avbp_sol,
                 ]
        put_files += [avbp_sol, avbp_mesh]
    else:
        raise XDRInterrupt("Unknown value for xor_flame_geo: " + flame_geo)

    n_tau_type = ds.getValue("xor_n_and_tau", "modify")
    if n_tau_type == "single_values":
        d = ds.getListDict("values", "xor_n_and_tau")
        n_tau_path = "tmp.dat"
        with open(n_tau_path, 'w') as f:
            f.write("{0} {1} {2}".format(
                d["Frequency"], d["N"], d["tau"]))
    else:
        n_tau_path = ds.getValue("n_tau_data_path")

    n_type = ds.getValue("xor_n_type", "modify")
    if n_type == "global":
        script += ["se nt", n_type, basename(n_tau_path)]
    elif n_type == "crocco":
        area = ds.getValue("area", "xor_n_type", "modify")
        p_mean = ds.getValue("p_mean", "xor_n_type", "modify")
        gamma = ds.getValue("gamma", "xor_n_type", "modify")
        script += ["se nt", n_type, basename(n_tau_path), "{0} {1} {2}".format(area, p_mean, gamma)]
    elif n_type == "adim":
        q_bar = ds.getValue("q_bar", "xor_n_type", "modify")
        u_bar = ds.getValue("u_bar", "xor_n_type", "modify")
        script += ["se nt", n_type, basename(n_tau_path), q_bar + " " + u_bar]
    script += ["se re",
               " ".join(ds.getListValue("ptref_list", "modify")[1::2]),
               " ".join(ds.getListValue("vecref_list", "modify")[1::2])]
    put_files += [n_tau_path]

    script += ["wr fu"]
    script += ["qu\n"]

    pr.execute_script("-create_flame-",
                      "\n".join(script), put_files=put_files, 
                      get_files=[pr.libobj_file(label)])

def update_flame_params(pr):
    """Update the 'modify' tab with values from a selected existing flame"""
    ds = pr.ds
    selected_flame = ds.getValue("cho_update_flame_params")
    metas = pr.get_metas(paths=[pr.libobj_file(selected_flame)])[selected_flame]
    # Rebuild n-tau
    ds.removeNode("n_tau")
    ds.addChild("n_tau", "", "modify")
    single_val = len(metas["n2_tau"][0]) == 1
    if single_val:
        ds.addChild("xor_n_and_tau", "single_values", "n_tau")
        ds.addChild("single_values", "", "xor_n_and_tau")
        ds.addChild("values",
                    to_ds_list("Frequency N tau".split(),
                               [x for y in metas["n2_tau"] for x in y]),
                    "single_values")
    else:
        n_tau_file = metas["name"] + ".n2_tau.dat"
        with open(n_tau_file, "w") as f:
            f.write("# Freq N2 tau\n")
            np.savetxt(f, np.array(metas["n2_tau"]).T)
        ds.addChild("xor_n_and_tau", "from_file", "n_tau")
        ds.addChild("n_tau_data_path", n_tau_file, "xor_n_and_tau")

    # Rebuild xor_ndim
    ds.removeNode("xor_ndim")
    ndim = "three_d" if "3D" in metas["generation_method"] else "two_d"
    ds.addChild("xor_ndim", ndim, "modify")
    ds.addChild(ndim, "", "xor_ndim")
    ds.addChild("ptref", "", ndim)
    ds.addChild("ptref_list", to_ds_list('xyz', metas["pt_ref"]), "ptref")
    ds.addChild("vecref_list", to_ds_list('xyz', metas["vec_ref"]), "ptref")
    ds.addChild("xor_flame_geo", metas["generation_method"], ndim)
    ds.addChild(metas["generation_method"], "", "xor_flame_geo")
    sh_args = json.loads(metas["shape_params"])[1]
    if metas["generation_method"] == "analytical2D_disc":
        ds.addChild("ana_flame_center", to_ds_list('xy', sh_args[0]), metas["generation_method"])
        ds.addChild("ana_flame_radius", str(sh_args[1]), metas["generation_method"])
    elif metas["generation_method"] == "analytical2D_rectangle":
        ds.addChild("ana_flame_pt_min", to_ds_list('ptmin_x ptmin_y'.split(), sh_args[0]), metas["generation_method"])
        ds.addChild("ana_flame_pt_max", to_ds_list('ptmax_x ptmax_y'.split(), sh_args[1]), metas["generation_method"])
    elif metas["generation_method"] == "analytical3D_brick":
        ds.addChild("ana_flame_pt_min", to_ds_list('ptmin_x ptmin_y ptmin_z'.split(), sh_args[0]), metas["generation_method"])
        ds.addChild("ana_flame_pt_max", to_ds_list('ptmax_x ptmax_y ptmax_z'.split(), sh_args[1]), metas["generation_method"])
    elif metas["generation_method"] == "analytical3D_sphere":
        ds.addChild("ana_flame_center", to_ds_list('xyz', sh_args[0]), metas["generation_method"])
        ds.addChild("ana_flame_radius", str(sh_args[1][0]), metas["generation_method"])
    elif metas["generation_method"] == "analytical3D_cylinder":
        ds.addChild("ana_flame_center", to_ds_list('xyz', sh_args[0]), metas["generation_method"])
        ds.addChild("ana_flame_radius", str(sh_args[1][0]), metas["generation_method"])
        ds.addChild("ana_flame_vector", to_ds_list('u_x u_y u_z'.split(), sh_args[2]), metas["generation_method"])
    elif metas["generation_method"] == "avbp_scalar_threshold_3D":
        ds.addChild("avbp_sol", metas["avbp_sol"])
        ds.addChild("avbp_mesh", metas["avbp_mesh"])
        ds.addChild("scal", metas["field"])
        ds.addChild("threshold", str(metas["threshold"]))
    else:
        raise ValueError("Unknown generation method {}. Cannot set default "
                         "values in 'modify' tab"
                         .format(metas["generation_method"]))

def to_ds_list(keys, values):
    out = []
    for k, v in zip(keys, values):
        out.extend([str(k), str(v)])
    return ";".join(out)

if __name__ == '__main__':
    import sys
    pr = FlameTransferProcess('dataset.xml')
    if "update" in sys.argv:
        pr.log.info("Updating interface, no modifications to the flame")
        update_flame_params(pr)
    else:
        process_modify(pr)
        pr.update_libobjs()
    pr.finish()

