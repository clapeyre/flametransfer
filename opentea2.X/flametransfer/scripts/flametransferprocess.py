#!/usr/bin/env python
"""
Concentrate all input/output methods for FlameTransfer scripts

Created June 16th, 2016 by C. Lapeyre (clapeyre@imft.fr, lapeyre@cerfacs.fr)
"""

import os
import logging
import json
import numpy as np

from os.path import join, isfile, isdir, basename

from XDR2 import LibProcess, RUN_CURRENT

__all__ = ['FlameTransferProcess']


class FlameTransferProcess(LibProcess):
    """Child class of LibProcess, adding FlameTransfer specific methods"""
    def __init__(self, name='dataset.xml'):
        LibProcess.__init__(self, name)
        self.log = logging.getLogger(__name__)
        self.debug = True

    def mesh_file(self, name): return "{}.mesh.h5".format(name)

    def libobj_file(self, name): return "{}.flame.h5".format(name)

    def libobj_name(self, path): return basename(path)[:-9]

    def update_libobjs(self):
        """Override mother method to update list_non_transformed_libobjs"""
        LibProcess.update_libobjs(self)
        metas = self.get_metas()
        non_transformed = [name for name in self.libobjs
                           if "transform" not in metas[name].keys()]
        self.ds.setValue(";".join(non_transformed),
                         "list_non_transformed_libobjs")

    def write_metas(self):
        metas = self.get_metas()
        for order, flame_name in enumerate(set(self.libobjs)
                                           - set(self.unwritten_libobjs)):
            ipcode = "item_{0:07d}".format(order)
            self.ds.addChild(ipcode, flame_name, "mul_libobjs")
            self.ds.addChild("written", "yes", ipcode, "mul_libobjs")
            self.ds.addChild(
                    "flame_details", metas[flame_name]["generation_method"],
                    ipcode, "mul_libobjs")
            details = "\n".join(("Freq, N2, tau",
                                 str(metas[flame_name]["n2_tau"]),
                                 ""
                                 "Ref point",
                                 str(metas[flame_name]["pt_ref"]),
                                 "Ref vector",
                                 str(metas[flame_name]["vec_ref"])))
            self.ds.addChild("ftf_details", details, ipcode, "mul_libobjs")
        for order, flame_name in enumerate(self.unwritten_libobjs):
            ipcode = "item_{0:07d}".format(order + len(self.libobjs))
            self.ds.addChild(ipcode, flame_name, "mul_libobjs")
            self.ds.addChild("written", "no", ipcode, "mul_libobjs")
            self.ds.addChild("flame_details", "not set", ipcode, "mul_libobjs")
            self.ds.addChild("ftf_details", "not set", ipcode, "mul_libobjs")

    def _rename(self, src_name, dest_name):
        """Specific implementation of flame renaming"""
        script = ["read " + self.libobj_file(src_name),
                  "set static name string",
                  dest_name,
                  "write flame"
                  ]
        self.execute_script(
                '-rename-', "\n".join(script),
                put_files=[self.libobj_file(src_name)],
                get_files=[self.libobj_file(dest_name)])

    def import_libobj(self, path):
        """Specific implementation of flame importing"""
        keys, values = self.get_metas(paths=[path]).items()[0]
        if values["name"] in self.libobjs + self.unwritten_libobjs:
            self.log.error(
                    values["name"] + " is already a declared flame. Rename "
                    "or delete it first.")
            raise ValueError
        self.copy_file(path, '.')
        self.add_libobj(values["name"])

    def duplicate_libobj(self, name, new_name=None):
        """Specific implementation of object duplication"""
        if new_name is None:
            new_name = name + '_duplicate'
        script = ["read " + self.libobj_file(name),
                  "set static name string",
                  new_name,
                  "wr fl"
                  ]
        self.execute_script(
                '-replicate-', "\n".join(script),
                put_files=[self.libobj_file(name)],
                get_files=[self.libobj_file(new_name)])
        self.add_libobj(new_name)

    def replicate_translate(self, name, vector, final_nb):
        """Replicate the flame by translating
        Translate by vector to obtain final_nb flames
        """
        script = ["read " + self.libobj_file(name)]
        new_flames = ["{0}_{1}".format(name, i) for i in range(2, final_nb+1)]
        translate = vector
        for target in new_flames:
            script += ["copy {0} {1}".format(name, target),
                       "transform translate",
                       " ".join(str(e) for e in translate),
                       "set static transform string",
                       "translate {0} by {1}".format(name, translate)]
            translate += vector
        script += ["wr fl"]
        self.execute_script(
                '-replicate-', "\n".join(script),
                put_files=[self.libobj_file(name)],
                get_files=[self.libobj_file(f) for f in new_flames])
        self.add_libobj(*new_flames)

    def replicate_rotate(self, name, angle, final_nb):
        """Replicate the flame by rotating
        Rotate by angle (degrees) to obtain final_nb flames
        """
        angles = [i*angle for i in range(1, final_nb)]
        script = ["read " + self.libobj_file(name)]
        new_flames = ["{0}_{1}".format(name, angl) for angl in angles]
        for angl, target in zip(angles, new_flames):
            script += ["copy {0} {1}".format(name, target),
                       "transform rotate",
                       "x",
                       str(angl),
                       "set static transform string",
                       "rotate {0} by {1}".format(name, angl)]
        script += ["wr fl"]
        self.execute_script(
                '-replicate-', "\n".join(script),
                put_files=[self.libobj_file(name)],
                get_files=[self.libobj_file(f) for f in new_flames])
        self.add_libobj(*new_flames)

    def execute_script(self, action, script, put_files=[], get_files=[]):
        extor = self.get_executor(action)
        with open(join(extor.execute_dir, 'script.ft'), 'w') as f:
            f.write(script)
            f.write('\nqu\n')
        self.log.debug("Content of script.ft:")
        for line in script.split('\n'):
            self.log.debug(line)
        for f in put_files:
            self.copy_file(f, extor.execute_dir)
        try:
            extor.execute()
        except RuntimeError, err:
            self.log.error("Flametransfer execution went wrong. "
                           "Check flametransfer.log for details")
            raise RuntimeError(err)
        for f in get_files:
            self.copy_file(join(extor.execute_dir, f), '.')
        extor.finish()

    def get_metas(self, paths=None):
        """Parse output of 'list metas'"""
        paths = paths if paths else [self.libobj_file(name)
                                     for name in self.libobjs]
        self.log.debug("Fetching metas for " + " ".join(paths))
        names = [self.libobj_name(path) for path in paths]
        script = ["read {} metas_only\nlist json".format(basename(f))
                  for f in paths]
        metas = ["{}.metas.json".format(self.libobj_name(path))
                 for path in paths]
        self.execute_script("-get_metas-",
                            "\n".join(script),
                            put_files=paths,
                            get_files=metas)
        dicts = {name: json.load(open(path))
                 for name, path in zip(names, metas)}
        [d.update({k: np.array(v) for k, v in d.items()
                   if k in d["numpies"]})
         for d in dicts.values()]
        if not self.debug:
            [os.remove(meta_file) for meta_file in metas]
        return dicts
