#!/usr/bin/env python
"""
Concentrate all input/output methods for FlameTransfer scripts

Created June 16th, 2016 by C. Lapeyre (clapeyre@imft.fr, lapeyre@cerfacs.fr)
"""

import logging
import os

from os.path import join, isfile, isdir, basename
from glob import glob

from XDR2 import LibProcess, RUN_CURRENT, COMMON
from XDR2.exceptions import XDRException

__all__ = ['FlameTransferProcess']

class FlameTransferProcess(LibProcess):
    """Child class of LibProcess, adding FlameTransfer specific methods"""
    def __init__(self, name='dataset.xml'):
        LibProcess.__init__(self, name)
        self.log = logging.getLogger(__name__)
        self.debug = True
        self.update_flames()

    def update_flames(self):
        """Update flame library multiple
        
        All .h5 files present should be accounted for in mul_flames.
        Any flame in mul_flames not present as an .h5 file is not yet written.
        It must stay in the multiple, but with details empty.
        """
        flame_files = glob(self.flame_file('*'))
        self.log.debug('Found h5 files of flames: ' + repr(flame_files))
        self.flames = [self.flame_name(path) for path in flame_files]
        self.ds.setValue(";".join(self.flames), "list_written_flames")
        for flame in self.flames:
            self.ds.addToUniqList(flame, "list_flames_in_project")
        self.unwritten_flames = [
                f for f in self.ds.getListValue("list_flames_in_project")
                if f not in self.ds.getListValue("list_written_flames")]
        self.log.debug("Unwritten flames in xml: " + repr(self.unwritten_flames))
        rebuild = (set(self.flames + self.unwritten_flames) != 
                   set([self.ds.getValue(node) for node in self.ds.getChildrenName("mul_flames")]))
        if rebuild:
            self.log.debug("Rebuilding multiple flame library")
            self.ds.removeNode("mul_flames")
            self.ds.addChild("mul_flames", "", "flame_library")
            metas = self.get_metas()
            for order, flame_name in enumerate(self.flames):
                ipcode="item_{0:07d}".format(order)
                self.ds.addChild(ipcode, flame_name, "mul_flames")
                self.ds.addChild(
                        "flame_details", metas[flame_name]["generation_method"],
                        ipcode, "mul_flames")
                self.ds.addChild(
                        "ftf_details",
                        "Freq, N2, tau\n" + metas[flame_name]["n2_tau"],
                        ipcode, "mul_flames")
            for order, flame_name in enumerate(self.unwritten_flames):
                ipcode="item_{0:07d}".format(order + len(self.flames))
                self.ds.addChild(ipcode, flame_name, "mul_flames")
                self.ds.addChild("flame_details", "not set", ipcode, "mul_flames")
                self.ds.addChild("ftf_details", "not set", ipcode, "mul_flames")

    def mesh_file(self, name): return "flame_{}.mesh.h5".format(name)

    def flame_file(self, name): return "flame_{}.sol.h5".format(name)

    def flame_name(self, path): return basename(path)[6:-7]

    def add_flame(self, name):
        if name in self.flames + self.unwritten_flames:
            self.log.error("Flame {} already exists".format(name))
            raise ValueError
        if isfile(self.flame_file(name)):
            self.flames.append(name)
            current_flame_list = self.ds.getValue("list_written_flames")
            self.ds.setValue(current_flame_list + ";" + name, "list_written_flames")
        else:
            self.unwritten_flames.append(name)
        current_flame_list = self.ds.getListValue("list_flames_in_project")
        self.ds.setValue(";".join(current_flame_list + [name]), "list_flames_in_project")
        self.update_flames()

    def rename_flame(self, src_name, dest_name):
        if src_name not in self.flames + self.unwritten_flames:
            self.log.error("No flame {} exists, cannot rename".format(src_name))
            raise XDRException
        elif src_name in self.flames:
            script = ["re " + self.flame_file(src_name),
                      "se me",
                      "name",
                      "st",
                      dest_name,
                      "wr fl"
                      ]
            self.execute_script(
                    '-rename-', "\n".join(script),
                    put_files=[self.flame_file(src_name)],
                    get_files=[self.flame_file(dest_name)])
        self.add_flame(dest_name)
        self.del_flame(src_name)

    def import_flame(self, path):
        keys, values = self.get_metas(paths=[path]).items()[0]
        if values["name"] in self.flames + self.unwritten_flames:
            self.log.error(
                    values["name"] + " is already a declared flame. Rename "
                    "or delete it first.")
            raise ValueError
        self.copy_file(path, '.')
        self.add_flame(values["name"])

    def duplicate_flame(self, name):
        script = ["re " + self.flame_file(name),
                  "se me",
                  "name",
                  "st",
                  name + '_duplicate',
                  "wr fl"
                  ]
        self.execute_script(
                '-duplicate-', "\n".join(script),
                put_files=[self.flame_file(name)],
                get_files=[self.flame_file(name + '_duplicate')])
        self.add_flame(name + '_duplicate')


    def del_flame(self, name):
        if name in self.flames:
            self.log.warning("Deleting flame file " + self.flame_file(name))
            self.flames.remove(name)
            os.remove(self.flame_file(name))
            flame_list = self.ds.getListValue("list_written_flames")
            flame_list.remove(name)
            self.ds.setValue(";".join(flame_list), "list_written_flames")
        elif name in self.unwritten_flames:
            self.log.info("Forgetting unwritten flame " + name)
        else:
            self.log.warning(name + " is not a known flame. Nothing to delete.")
            return
        flame_list = self.ds.getListValue("list_flames_in_project")
        flame_list.remove(name)
        self.ds.setValue(";".join(flame_list), "list_flames_in_project")
        for child in sorted(self.ds.getChildrenName("mul_flames")):
            if self.ds.getValue(child) == name:
                self.ds.removeNode(child)

    def execute_script(self, action, script, put_files=[], get_files=[]):
        extor = self.get_executor(action)
        with open(join(extor.execute_dir, 'script.ft'), 'w') as f:
            f.write(script)
            f.write('\nqu\n')
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
        """Parse output of 'list metas' """
        script = []
        metas = []
        paths = paths if paths else [self.flame_file(name) for name in self.flames]
        self.log.debug("Fetching metas for " + " ".join(paths))
        names = [self.flame_name(path) for path in paths]
        iterable = zip(names, paths)
        for name, flame_file in iterable:
            script += ["re {} metas_only".format(basename(flame_file)),
                       "li me metas_{}.dat".format(name)]
            metas.append("metas_{}.dat".format(name))
        script += ["qu"]
        self.execute_script("-get_metas-",
                            "\n".join(script),
                            put_files=paths,
                            get_files=metas)
        dicts = {}
        for meta_file, name in zip(metas, names):
            self.log.debug("Meta file " + meta_file)
            for line in open(meta_file, 'r').readlines():
                self.log.debug(line.strip())
            keys, values = zip(*[line.split('|') for line in open(meta_file, 'r').readlines()[2:]])
            keys = [k.strip() for k in keys]
            def value_handler(val):
                if '[' in val: val = val.translate(None, '[]\n')
                return val.strip()
            values = map(value_handler, values)
            if not self.debug: os.remove(meta_file)
            dicts[name] = dict(zip(keys, values))
        for meta_file in metas: os.remove(meta_file)
        return dicts
