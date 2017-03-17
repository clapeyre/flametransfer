#!/usr/bin/env python
"""
Database of flames and subsequent actions

Created Jan 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""

import os
import shutil

from glob import glob
from UserList import UserList
from textwrap import dedent

import numpy as np

from h5py import File

from activeflame import ActiveFlame
from geometry import Point, NormalVector
from constants import DEBUG
from tools import visu

class FlameDB(UserList, object):
    """List-like object to contain all flames

    Enables list-like access, plus supervision of actions performed on
    multiple flames.
    """
    def __init__(self, hip_wrapper):
        UserList.__init__(self)
        self.hip_wrapper = hip_wrapper
        self.current = None

    @property
    def names(self):
        return [fla.name for fla in self.data]

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

    def _get_mesh_voln(self, mesh):
        """Get volume at nodes.
        
        Read and write with modern hip to recompute if needed
        """
        fil = File(mesh, 'r')
        if "VertexData" not in fil.keys():
            fil.close()
            script = dedent("""
                se ch 0
                re hd -a {0}
                wr hd {1}
                qu""".format(mesh, "dummy"))
            self.hip_wrapper.execute(script)
            fil = File("dummy.mesh.h5", 'r')
        out = fil["VertexData/volume"].value
        fil.close()
        return out

    def import_avsp5x(self, mesh, sol):
        """Read all flames from target AVSP 5.X mesh and sol"""
        voln = self._get_mesh_voln(mesh)
        with File(sol, 'r') as fil:
            nnode = fil["Parameters/nnode"].value[0]
            zones = []
            refs = []
            for name, group in fil["/Zones"].items():
                mask = np.zeros((nnode,), dtype=bool)
                mask[group["znode->node"].value - 1] = True
                zones.append(mask)
                refs.append([group["SolParameters/Ptref"].value,
                             group["SolParameters/Unref"].value])
            n_local = fil["/Additionals/nu_global"].value
            integral = n_local*voln
            flame_flag = n_local > 1e-10
            tau = fil["/Additionals/tauu_global"].value
            flame_flags = [flame_flag*zone for zone in zones]
            taus = [tau[flag].mean() for flag in flame_flags]
            n_globals = [np.sum(integral[zone]) for zone in zones]
        for i in range(len(zones)):
            flame = ActiveFlame('tmp', self.hip_wrapper)
            flame.metas.name = sol.replace(".sol.h5", '') + "_{:03}".format(i + 1)
            flame.define_flame_scatter(mesh, flame_flags[i])
            flame.set_n2_tau(1., n_globals[i], taus[i])
            pt_ref, vec_ref = refs[i]
            flame.metas.pt_ref = Point(pt_ref)
            flame.metas.vec_ref = NormalVector(vec_ref)
            self.append(flame)
        self.current = self.data[-1]

    def import_avsp6x(self, sol):
        """Read all flames from target AVSP 6.X sol"""
        with File(sol, 'r') as fil:
            for name, group in fil["/Flames"].items():
                flame = ActiveFlame('tmp', self.hip_wrapper)
                flame._read_metas(group, number=int(name))
                self.append(flame)
        self.current = self.data[-1]

