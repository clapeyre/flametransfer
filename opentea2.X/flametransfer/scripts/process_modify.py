#!/usr/bin/env python
"""
process_modify.py

Process 'Modify' tab for FlameTransfer app

Created December 2016 by C. Lapeyre (lapeyre@cerfacs.fr)
"""

import json

from os.path import basename

import numpy as np

from opentea import OTInterrupt
from flametransferprocess import FlameTransferProcess


def process_modify(pro):
    """Modify existing or new virgin flame"""
    dst = pro.ds

    label = dst.getValue("cho_write_to_flame")
    put_files = []
    flame_geo = dst.getValue("xor_flame_geo", "modify")
    if flame_geo == "analytical2D_disc":
        center = dst.getListValue("ana_flame_center", "modify")[1::2]
        radius = dst.getValue("ana_flame_radius", "modify")
        script = ["ge disc " + label,
                  "{0[0]} {0[1]}".format(center),
                  "{}".format(radius)]
    elif flame_geo == "analytical2D_rectangle":
        pt_min = dst.getListValue("ana_flame_pt_min", "modify")[1::2]
        pt_max = dst.getListValue("ana_flame_pt_max", "modify")[1::2]
        script = ["ge rect " + label]
        script += ["{0[0]} {0[1]}".format(vect)
                   for vect in [pt_min, pt_max]]
    elif flame_geo == "analytical3D_brick":
        pt_min = dst.getListValue("ana_flame_pt_min", "modify")[1::2]
        pt_max = dst.getListValue("ana_flame_pt_max", "modify")[1::2]
        script = ["ge brick " + label]
        script += ["{0[0]} {0[1]} {0[2]}".format(vect)
                   for vect in [pt_min, pt_max]]
    elif flame_geo == "analytical3D_sphere":
        center = dst.getListValue("ana_flame_center", "modify")[1::2]
        radius = dst.getValue("ana_flame_radius", "modify")
        script = ["ge sphere " + label,
                  "{0[0]} {0[1]} {0[2]}".format(center),
                  "{}".format(radius)]
    elif flame_geo == "analytical3D_cylinder":
        center = dst.getListValue("ana_flame_center", "modify")[1::2]
        radius = dst.getValue("ana_flame_radius", "modify")
        vector = dst.getListValue("ana_flame_vector", "modify")[1::2]
        pro.log.debug("Creating cylinder with center " + repr(center))
        pro.log.debug("                       radius " + repr(radius))
        pro.log.debug("                       vector " + repr(vector))
        script = ["ge cylinder " + label,
                  "{0[0]} {0[1]} {0[2]}".format(center),
                  "{}".format(radius),
                  "{0[0]} {0[1]} {0[2]}".format(vector)]
    elif flame_geo in ["avbp_scalar_threshold_3D", "avbp_scalar_threshold_2D"]:
        avbp_sol = dst.getValue("avbp_sol", "xor_flame_geo", "modify")
        avbp_mesh = dst.getValue("avbp_mesh", "xor_flame_geo", "modify")
        scal = dst.getValue("scal", "xor_flame_geo", "modify")
        thresh = dst.getValue("thresh", "xor_flame_geo", "modify")
        script = ["ge avbp " + label,
                  basename(avbp_mesh),
                  basename(avbp_sol),
                  scal,
                  thresh,
                  "set static avbp_mesh string",
                  avbp_mesh,
                  "set static avbp_sol string",
                  avbp_sol]
        put_files += [avbp_sol, avbp_mesh]
    else:
        raise OTInterrupt("Unknown value for xor_flame_geo: " + flame_geo)

    n_tau_type = dst.getValue("xor_n_and_tau", "modify")
    if n_tau_type == "single_values":
        d = dst.getListDict("values", "xor_n_and_tau")
        n_tau_path = "tmp.dat"
        with open(n_tau_path, 'w') as f:
            f.write("{0} {1} {2}".format(
                d["Frequency"], d["N"], d["tau"]))
    else:
        n_tau_path = dst.getValue("n_tau_data_path")

    n_type = dst.getValue("xor_n_type", "modify")
    if n_type == "global":
        script += ["se nt", n_type, basename(n_tau_path)]
    elif n_type == "crocco":
        area = dst.getValue("area", "xor_n_type", "modify")
        p_mean = dst.getValue("p_mean", "xor_n_type", "modify")
        gamma = dst.getValue("gamma", "xor_n_type", "modify")
        script += ["se nt",
                   n_type,
                   basename(n_tau_path),
                   "{0} {1} {2}".format(area, p_mean, gamma)]
    elif n_type == "adim":
        q_bar = dst.getValue("q_bar", "xor_n_type", "modify")
        u_bar = dst.getValue("u_bar", "xor_n_type", "modify")
        script += ["se nt", n_type, basename(n_tau_path), q_bar + " " + u_bar]
    script += ["se re",
               " ".join(dst.getListValue("ptref_list", "modify")[1::2]),
               " ".join(dst.getListValue("vecref_list", "modify")[1::2])]
    put_files += [n_tau_path]

    script += ["wr fu"]
    script += ["qu\n"]

    pro.execute_script("-create_flame-",
                       "\n".join(script), put_files=put_files,
                       get_files=[pro.libobj_file(label)])


def update_flame_params(pro):
    """Update the 'modify' tab with values from a selected existing flame"""
    dst = pro.ds
    selected_flame = dst.getValue("cho_update_flame_params")
    metas = pro.get_metas(
        paths=[pro.libobj_file(selected_flame)])[selected_flame]
    gen_meth = metas["generation_method"]
    # Rebuild n-tau
    dst.removeNode("n_tau")
    dst.addChild("n_tau", "", "modify")
    single_val = len(metas["n2_tau"][0]) == 1
    if single_val:
        dst.addChild("xor_n_and_tau", "single_values", "n_tau")
        dst.addChild("single_values", "", "xor_n_and_tau")
        dst.addChild("values",
                     to_ds_list("Frequency N tau".split(),
                                [x for y in metas["n2_tau"] for x in y]),
                     "single_values")
    else:
        n_tau_file = metas["name"] + ".n2_tau.dat"
        with open(n_tau_file, "w") as f:
            f.write("# Freq N2 tau\n")
            np.savetxt(f, np.array(metas["n2_tau"]).T)
        dst.addChild("xor_n_and_tau", "from_file", "n_tau")
        dst.addChild("n_tau_data_path", n_tau_file, "xor_n_and_tau")

    # Rebuild geometry
    dst.removeNode("three_d")
    dst.addChild("three_d", "", "modify")
    dst.addChild("xor_flame_geo", "", "three_d")
    dst.addChild(gen_meth, "", "xor_flame_geo")
    sh_args = json.loads(metas["shape_params"])[1]
    print sh_args
    if gen_meth == "analytical3D_brick":
        dst.addChild("ana_flame_pt_min",
                     to_ds_list('ptmin_x ptmin_y ptmin_z'.split(), sh_args[0]),
                     gen_meth)
        dst.addChild("ana_flame_pt_max",
                     to_ds_list('ptmax_x ptmax_y ptmax_z'.split(), sh_args[1]),
                     gen_meth)
    elif gen_meth == "analytical3D_sphere":
        dst.addChild("ana_flame_center",
                     to_ds_list('xyz', sh_args[0]),
                     gen_meth)
        dst.addChild("ana_flame_radius", str(sh_args[1][0]), gen_meth)
    elif gen_meth == "analytical3D_cylinder":
        dst.addChild("ana_flame_center",
                     to_ds_list('xyz', sh_args[0]),
                     gen_meth)
        dst.addChild("ana_flame_radius", str(sh_args[1][0]), gen_meth)
        dst.addChild("ana_flame_vector",
                     to_ds_list('u_x u_y u_z'.split(), sh_args[2]),
                     gen_meth)
    elif gen_meth == "avbp_scalar_threshold_3D":
        dst.addChild("avbp_sol", metas["avbp_sol"], gen_meth)
        dst.addChild("avbp_mesh", metas["avbp_mesh"], gen_meth)
        dst.addChild("scal", metas["field"], gen_meth)
        dst.addChild("threshold", str(metas["threshold"]), gen_meth)
    else:
        raise ValueError("Unknown generation method {}. Cannot set default "
                         "values in 'modify' tab"
                         .format(gen_meth))


def to_ds_list(keys, values):
    out = []
    for k, v in zip(keys, values):
        out.extend([str(k), str(v)])
    return ";".join(out)


if __name__ == '__main__':
    import sys
    pro = FlameTransferProcess('dataset.xml')
    if "update" in sys.argv:
        pro.log.info("Updating interface, no modifications to the flame")
        update_flame_params(pro)
    else:
        process_modify(pro)
        pro.update_libobjs()
    pro.finish()
