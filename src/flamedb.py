#!/usr/bin/env python
"""
Database of flames and subsequent actions

Created Jan 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""

import os
import shutil

from UserList import UserList

import numpy as np

from h5py import File

from activeflame import ActiveFlame

class FlameDB(UserList, object):
    current = None

    def export(self, mesh, sol, target):
        """Export all flames to target"""
        if target == "avsp5.6":
            shutil.copy(sol, "avsp_tmp0.sol.h5")
            with File("avsp_tmp0.sol.h5", 'a') as f:
                if not "Average" in f.keys():
                    f.create_group("Average")
                    f["Average/flames"] = np.zeros(f["Parameters/nnode"].value)
                elif not "flames" in f["Average"].keys():
                    f["Average/flames"] = np.zeros(f["Parameters/nnode"].value)
            for i,flame in enumerate(self.data):
                print " --- Exporting flame {0} with number {1}".format(flame.name, i+1)
                flame.export_avsp(mesh, "avsp_tmp{}.sol.h5".format(i), i+1)
                os.rename("avsp.sol.h5", "avsp_tmp{}.sol.h5".format(i+1))
            os.rename("avsp_tmp{}.sol.h5".format(i+1), "avsp.sol.h5")
            with File("avsp.sol.h5", 'a') as avsp:
                for i,flame in enumerate(self.data):
                    flame.write_group(avsp, number=i+1)
        else:
            raise AssertionError("Unknown target " + target)

    def import_flames(self, sol, target, hip_exec):
        """Read all flames from target"""
        if target == "avsp5.6":
            with File(sol, 'r') as f:
                for name, group in f["/Flames"].items():
                    flame = ActiveFlame('tmp', hip_exec)
                    flame._read_metas(group, number=int(name))
                    self.append(flame)
            self.current = self.data[-1]
        else:
            raise AssertionError("Unknown target " + target)

