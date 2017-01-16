#!/usr/bin/env python
"""
ActiveFlame class to handle all information pertaining to an active flame.
Handle Read / Write of H5 Flame files.

Created November 2016 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""
DEBUG = False

import os
import sys
import time
import logging
import subprocess
import commands
import copy
import json
import numpy as np

from os.path import join, isdir
from math import pi
from h5py import File
from textwrap import dedent
from StringIO import StringIO
from tempfile import TemporaryFile

from geometry import (Rectangle, Sphere, Cylinder, Brick, Disc,
                      ScatterShape2D, ScatterShape3D, NormalVector, Vector,
                      Point, json2shape, shape2json)
from flamemetas import FlameMetas
from constants import VERSION

class ActiveFlame(object):
    """Flame holder class associated to a flame hdf5 file"""
    def __init__(self, name, hip_exec, **kwargs):
        self.log = logging.getLogger(__name__)
        self.debug = True
        self.meshpoints = None
        self.shape = None
        self.metas = FlameMetas(name=name)
        self.hip_exec = hip_exec
        self.last_hip_output = ""
        self.__dict__.update(kwargs)

    @property
    def mesh_file(self):
        return 'flame_{0}.mesh.h5'.format(self.metas.name)

    @property
    def flame_file(self):
        return 'flame_{0}.sol.h5'.format(self.metas.name)


    def __getstate__(self):
        d = dict(self.__dict__)
        del d['log']
        return d

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.log = logging.getLogger(__name__)

    def exec_hip(self, script):
        """Execute the current hip script."""
        script += "\nqu\n"
        path = 'script.hip'
        self.last_hip_output = ""
        self.log.debug("Executing hip script:") 
        [self.log.debug(" > " + line) for line in script.split('\n')]
        with open(path, 'w') as f:
            f.write(script)
        with TemporaryFile() as output:
            process = subprocess.Popen(
                    [self.hip_exec, path],
                    stdin=subprocess.PIPE,
                    stdout=output,
                    stderr=subprocess.STDOUT
                    )
            while process.poll() is None:
                where = output.tell()
                lines = output.read()
                if not lines:
                    time.sleep(0.1)
                    output.seek(where)
                else:
                    self.last_hip_output += lines
            self.last_hip_output += output.read()
            if process.wait() != 0:
                print "*** error in hip script. See log"
            self.log.debug("Hip execution log")
            for line in self.last_hip_output.split('\n'):
                self.log.debug(line)
        print "\n --- Done executing hip"

    def compute_n_crocco(self, n1, area, p_mean, gamma):
        """Define gain (N2) using Crocco's analytical formulation"""
        return gamma * p_mean * area / (gamma - 1.) * n1

    def set_n1_tau(self, freq, tau, n1, area, p_mean, gamma):
        """Set N and tau from Crocco model (N1)"""
        self.set_n2_tau(freq, tau, self.compute_n_crocco(n1, area, p_mean, gamma))

    def set_n2_tau(self, freq, tau, n2):
        """Set N and tau from global values (same as AVSP internal)"""
        self.metas.n2_tau = np.array((freq, n2, tau), ndmin=2)
        # Ugly!! There's probably a better way to keep orientation consistent
        if self.metas.n2_tau.shape[0] != 3:    
            self.metas.n2_tau = self.metas.n2_tau.T

    def set_n3_tau(self, freq, tau, n3, u_bar, q_bar):
        """Set N and tau from non dimensional values (N3)"""
        self.u_bar = u_bar
        self.q_bar = q_bar
        self.set_n2_tau(freq, tau, n3 * u_bar / q_bar)

    def get_n3(self):
        """Get value of N3 using u_bar and q_bar"""
        if hasattr(self.metas, "q_bar"):
            return self.n2 * self.metas.q_bar / self.metas.u_bar
        else:
            self.log.error("q_bar info not available, cannot compute n3")

    def define_threshold_flame(self, avbp_mesh, avbp_sol, field, thresh):
        """Define flame from avbp scalar and threshold"""
        with File(avbp_mesh, 'r') as mesh, File(avbp_sol, 'r') as sol:
            def find_field(name):
                if field.lower() in name.lower(): return name
            selected_field = sol.visit(find_field)
            if selected_field is None:
                print "*** no field containing {0} found in {1}".format(field, avbp_sol)
                return
            data = sol[sol.visit(find_field)].value
            above = data > np.float(thresh)
            x = mesh["Coordinates/x"].value
            y = mesh["Coordinates/y"].value
            x_in = x[above]
            y_in = y[above]
            pt_min = [x_in.min(), y_in.min()]
            pt_max = [x_in.max(), y_in.max()]
            klass = ScatterShape2D
            self.metas.generation_method = "avbp_scalar_threshold_2D"
            if 'z' in mesh["Coordinates"].keys():
                self.metas.generation_method = "avbp_scalar_threshold_3D"
                z = mesh["Coordinates/z"].value
                z_in = z[above]
                pt_min += [z_in.min()]
                pt_max += [z_in.max()]
                klass = ScatterShape3D
        self.log.debug("Generating scattershape from AVBP")
        self._shape_metas(klass(pt_min, pt_max))
        self.make_mesh()
        with File("dummy_avbp.h5", 'w') as f:
            f.create_group
            f.create_group("/Additionals")
            f['/Additionals/n_tau_flag'] = 1.0*above
            f.create_group("/Parameters")
            f['/Parameters/ndim'] = np.array([self.metas.ndim])
            f['/Parameters/nnode'] = np.array([self.meshpoints.size])
            f['/Parameters/versionstring'] = "flametransfer_v" + self.metas.version
        with File("dummy_flametrans.h5", 'w') as f:
            f.create_group
            f.create_group("/Additionals")
            f['/Additionals/n_tau_flag'] = 0.0*above
            f.create_group("/Parameters")
            f['/Parameters/ndim'] = np.array([self.metas.ndim])
            f['/Parameters/nnode'] = np.array([self.meshpoints.size])
            f['/Parameters/versionstring'] = "flametransfer_v" + self.metas.version
        script = [
                "re hd -a {0} -s dummy_avbp.h5".format(avbp_mesh),
                "re hd -a {0} -s dummy_flametrans.h5".format(self.mesh_file),
                "set in-rim 0.1",
                "in gr 1",
                "wr hd dummy_flame",
                "qu",
                ]
        self.exec_hip('\n'.join(script))
        self._read_scatter_pts("dummy_flame.sol.h5")
        if not self.debug:
            os.remove("dummy_avbp.h5")
            os.remove("dummy_flametrans.h5")
            os.remove("dummy_flame.sol.h5")
            os.remove("dummy_flame.mesh.h5")
            os.remove("dummy_flame.asciiBound")

    def _read_scatter_pts(self, path):
        """Read scatter points from flame file"""
        with File(path, 'r') as f:
            self.inside_points = f['Additionals/n_tau_flag'].value

    def define_flame_disc(self, center, radius):
        """Define flame as disc"""
        self.metas.generation_method = "analytical2D_disc"
        self._make_analytic_shape(Disc(center, radius))

    def define_flame_rectangle(self, pt_min, pt_max):
        """Define flame as rectangle"""
        self.metas.generation_method = "analytical2D_rectangle"
        self._make_analytic_shape(Rectangle(pt_min, pt_max))

    def define_flame_sphere(self, center, radius):
        """Define spherical flame"""
        self.metas.generation_method = "analytical3D_sphere"
        self._make_analytic_shape(Sphere(center, radius))

    def define_flame_cylinder(self, center, radius, vector):
        """Define cylindrical flame"""
        self.log.debug("Defining cylinder flame")
        self.metas.generation_method = "analytical3D_cylinder"
        self._make_analytic_shape(Cylinder(center, radius, vector))

    def define_flame_brick(self, pt_min, pt_max):
        """Define flame as brick"""
        self.metas.generation_method = "analytical3D_brick"
        self._make_analytic_shape(Brick(pt_min, pt_max))

    def _shape_metas(self, shape):
        self.metas.ndim = shape.ndim
        self.metas.pt_min = copy.deepcopy(shape.vects["pt_min"])
        self.metas.pt_max = copy.deepcopy(shape.vects["pt_max"])
        self.metas.volume = shape.volume
        self.metas.shape_params = shape2json(shape)

    def _make_analytic_shape(self, shape):
        self._shape_metas(shape)
        self.make_mesh()
        self.inside_points = 1.0*shape.is_inside(self.meshpoints)

    def make_mesh(self):
        """Create and read mesh for flame
        
        Use bounding box to create box mesh containing flame geometry
        Then read resulting meshoints"""
        self.exec_hip(self._get_hip_script_generate())
        self._read_meshpoints()

    def _read_meshpoints(self):
        """Read ouput hdf5 meshfile and store coordinates"""
        with File(self.mesh_file, 'r') as f:
            if self.metas.ndim == 2:
                self.meshpoints = np.array((f['/Coordinates/x'].value,
                                            f['/Coordinates/y'].value)).T
            elif self.metas.ndim == 3:
                self.meshpoints = np.array((f['/Coordinates/x'].value,
                                            f['/Coordinates/y'].value,
                                            f['/Coordinates/z'].value)).T
            else:
                raise ValueError("ndim must be defined")
            self.log.debug("Mesh data stored with {0} points".format(self.meshpoints[:,0].size))

    def _get_hip_script_generate(self):
        """Write hip script for mesh generation"""
        assert self.metas.pt_min is not None, "Bounding box needed for mesh gen"
        line_3d = (
                "co 3d {0[2]} {1[2]} {2} z".format(self.metas.pt_min,
                                                   self.metas.pt_max,
                                                   self.metas.grid_size - 1)
                if self.metas.ndim == 3 else "")
        hip_script = dedent("""\
          ge {0[0]} {0[1]} {1[0]} {1[1]} {2} {2}
          {3}
          wr hd ./flame_{4}
          qu
          """).format(self.metas.pt_min, self.metas.pt_max,
                      self.metas.grid_size, line_3d, self.metas.name)
        return hip_script

    def transform(self, trans, *args):
        """Apply transformation to flame"""
        if trans == "t":
            self.metas.transform("translate", *args)
        elif trans == "s":
            self.metas.transform("scale", *args)
        elif trans == "r":
            self.metas.transform("rotate", *args)
            # Special case: box must be rotated, then interpolated back on a
            # cartesian grid aligned with x,y[,z]
            self._rotate(*args)
        else:
            raise ValueError
        self.make_mesh()

    def _rotate(self, axis, angle, degrees=True):
        """Write and rotate flame in hip, then recreate x,y[,z] aligned mesh"""
        self.write_h5(with_mesh=True)
        if not degrees: angle *= 180./np.pi
        rotate_script = dedent("""\
          re hd -a flame_{0}.mesh.h5 -s flame_{0}.sol.h5
          tr ro {1} {2}
          wr hd ./flame_{0}
          """.format(self.metas.name, axis, angle))
        self.exec_hip(rotate_script)

        self._read_scatter_pts(self.flame_file)
        self._read_meshpoints()
        flamepoints = self.meshpoints[self.inside_points > 0.01]
        self.metas.pt_min = Point(flamepoints.min(0))
        self.metas.pt_max = Point(flamepoints.max(0))

        interp_script = ["re hd -a flame_{0}.mesh.h5 -s flame_{0}.sol.h5".format(self.metas.name)]
        interp_script += self._get_hip_script_generate().split("\n")
        interp_script.insert(3, "in gr 1")
        self.exec_hip('\n'.join(interp_script))

    def export_avsp(self, avsp_mesh, avsp_sol):
        """Export flame to AVSP solution"""
        self.write_h5(with_mesh=True)
        self.exec_hip(self._get_hip_script_export_avsp(avsp_mesh, avsp_sol))

    def _get_hip_script_export_avsp(self, avsp_mesh, avsp_sol):
        """Generate hip script for flame interpolation on AVSP mesh"""
        hip_script = dedent("""\
          se ch 0
          se in-rim 0.2
          re hd -a flame_{0}.mesh.h5 -s flame_{0}.sol.h5
          re hd -a {1} -s {2}
          in gr 1
          zone add secteur 
          zone 1 element add all 
          zo 1 solparam add Ptref vec {3}
          zo 1 solparam add Unref vec {4}
          wr hd avsp_sol
          ex
          """).format(self.metas.name, avsp_mesh, avsp_sol,
                      " ".join(str(s) for s in self.metas.pt_ref),
                      " ".join(str(s) for s in self.metas.vec_ref),
                      )
        return hip_script

    def write_h5(self, with_metas=True, with_mesh=False):
        """Write h5 file for flame"""
        assert self.meshpoints is not None
        with File(self.flame_file, 'w') as f:
            f.create_group("/Additionals")
            if with_metas: self.metas.write_h5(f["/Additionals"].attrs)
            f['/Additionals/n_tau_flag'] = self.inside_points
            f.create_group("/Parameters")
            f['/Parameters/ndim'] = np.array([self.metas.ndim])
            f['/Parameters/nnode'] = np.array([self.meshpoints.size])
            f['/Parameters/versionstring'] = "flametransfer_v" + self.metas.version
            self.log.info("Wrote flame file to " + self.flame_file)
        if with_mesh:
            self.make_mesh()
            script = [
                    "re hd -a {0.mesh_file} -s {0.flame_file}".format(self),
                    "wr hd flame_{0}".format(self.metas.name),
                    "qu\n"
                    ]
            self.exec_hip("\n".join(script))
            if with_metas:
                with File(self.flame_file, 'a') as f:
                    self.metas.write_h5(f["/Additionals"].attrs)
            os.remove("flame_{0}.asciiBound".format(self.metas.name))
            self.log.info(
                    "Wrote flame and mesh files: flame_{0}.mesh.h5, "
                    "flame_{0}.mesh.xmf, flame_{0}.sol.h5".format(self.metas.name))

    def write_n_tau(self):
        """Write ascii n_tau file for flame"""
        assert self.metas.n2_tau is not None
        path = 'n_tau_{}.dat'.format(self.metas.name)
        buf = StringIO()
        np.savetxt(buf, self.metas.n2_tau.T)
        with open(path, 'w') as f:
            f.write("# Frequency N2 tau  - Generated using flametransfer\n")
            f.write(''.join([c if c.lower() != 'e' else 'd' for c in buf.getvalue()]))

    def load(self, path):
        """Load a flame at path"""
        self.load_metas(path)
        self.log.info("Loaded flame named " + self.metas.name)
        self.log.debug(repr(self.metas.__dict__))
        self.make_mesh()
        self._read_scatter_pts(path)

    def load_metas(self, path):
        """Load only the metas of flame at path"""
        self.log.debug("Loading metas for file " + path)
        with File(path, 'r') as flame:
            self.metas.load(flame['Additionals'].attrs)
        if float(self.metas.version) > float(VERSION):
            self.log.error("Cannot read future version flame. Please update your FlameTransfer app")
            raise ValueError

