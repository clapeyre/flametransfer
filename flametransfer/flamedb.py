#!/usr/bin/env python
"""
Database of flames and subsequent actions

Created Jan 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""

import os
import shutil
import logging

from glob import glob
from UserList import UserList
from textwrap import dedent

import numpy as np

from h5py import File

import geometry as geo
from activeflame import ActiveFlame
from constants import DEBUG
from tools import visu


class FlameDB(UserList, object):
    """List-like object to contain all flames

    Enables list-like access, plus supervision of actions performed on
    multiple flames.
    """
    def __init__(self, hip_wrapper):
        UserList.__init__(self)
        self.log = logging.getLogger(__name__)
        self.hip_wrapper = hip_wrapper
        self.current = None

    @property
    def names(self):
        """List of names of all the flames"""
        return [fla.name for fla in self.data]

    def export_avsp6x(self, mesh, sol):
        """Export all flames to target"""
        shutil.copy(sol, "avsp_tmp0.sol.h5")
        with File("avsp_tmp0.sol.h5", 'a') as fil:
            if "Average" not in fil.keys():
                fil.create_group("Average")
                fil["Average/flames"] = np.zeros(fil["Parameters/nnode"].value)
            elif "flames" not in fil["Average"].keys():
                fil["Average/flames"] = np.zeros(fil["Parameters/nnode"].value)
        for i, flame in enumerate(self.data):
            print (" --- Exporting flame {0} with number {1}"
                   .format(flame.name, i+1))
            flame.export_avsp(mesh, "avsp_tmp{}.sol.h5".format(i), i+1)
        shutil.copy("avsp_tmp{}.sol.h5".format(i+1), "avsp.sol.h5")
        visu(mesh, "avsp.sol.h5")
        with File("avsp.sol.h5", 'a') as avsp:
            for i, flame in enumerate(self.data):
                flame.write_group(avsp, number=i+1)
        if not DEBUG:
            for fil in glob("avsp_tmp*"):
                os.remove(fil)

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
            for _, group in fil["/Zones"].items():
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
            flame.metas.name = (sol.replace(".sol.h5", '')
                                + "_{:03}".format(i + 1))
            flame.define_flame_scatter(mesh, flame_flags[i])
            flame.set_n2_tau(1., n_globals[i], taus[i])
            pt_ref, vec_ref = refs[i]
            flame.metas.pt_ref = geo.Point(pt_ref)
            flame.metas.vec_ref = geo.NormalVector(vec_ref)
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

    def test(self):
        """Run simple test for FlameTransfer"""
        assert not self, "please run `test` at startup"
        msg = "tests failed. See flametransfer.log"
        try:
            cyl = ActiveFlame('test_flame', self.hip_wrapper)
            cyl.define_flame_cylinder([0, 0, 0], 1, [1, 1, 0])
            cyl.set_n2_tau([0, 100], [500, 1000], [0.001]*2)
            cyl.metas.pt_ref = geo.Point([0, 0, 0])
            cyl.metas.vec_ref = geo.NormalVector([1, 1, 0])
            self.append(cyl)
            self[-1].write_full()
        except Exception as err:
            self.log.error(err.message)
            raise AssertionError(msg)
        assert os.path.isfile("test_flame.flame.h5"), msg
        assert os.path.isfile("test_flame.flame.xmf"), msg
        assert os.path.isfile("test_flame.mesh.h5"), msg
        print " --- All tests passed."
