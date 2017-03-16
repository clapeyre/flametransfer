#!/usr/bin/env python
"""
Database of flames and subsequent actions

Created Jan 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""

import os
import shutil

from glob import glob
from UserList import UserList

import numpy as np

from h5py import File

from activeflame import ActiveFlame
from constants import DEBUG
from tools import visu

class FlameDB(UserList, object):
    """List-like object to contain all flames

    Enables list-like access, but also supervision of actions performed on
    all flames.
    """
    def __init__(self, hip_wrapper):
        UserList.__init__(self)
        self.hip_wrapper = hip_wrapper
        self.current = None

    def export_avsp6x(self, mesh, sol):
        """Export all flames to target"""
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
        shutil.copy("avsp_tmp{}.sol.h5".format(i+1), "avsp.sol.h5")
        visu(mesh, "avsp.sol.h5")
        with File("avsp.sol.h5", 'a') as avsp:
            for i,flame in enumerate(self.data):
                flame.write_group(avsp, number=i+1)
        if not DEBUG: [os.remove(fil) for fil in glob("avsp_tmp*")]

    def import_avsp5x(self, mesh, sol):
        """Read all flames from target AVSP 5.X target"""
        dummy = ActiveFlame('tmp', self.hip_wrapper)
        s
        with File(sol, 'r') as f:
            for name, group in f["/Flames"].items():
                flame = ActiveFlame('tmp', self.hip_wrapper)
                flame._read_metas(group, number=int(name))
                self.append(flame)
        self.current = self.data[-1]

    def import_avsp6x(self, sol):
        """Read all flames from target AVSP 6.X target"""
        with File(sol, 'r') as f:
            for name, group in f["/Flames"].items():
                flame = ActiveFlame('tmp', self.hip_wrapper)
                flame._read_metas(group, number=int(name))
                self.append(flame)
        self.current = self.data[-1]

