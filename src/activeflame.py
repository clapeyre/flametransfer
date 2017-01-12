#!/usr/bin/env python
"""
ActiveFlame class to handle all information pertaining to an active flame.
Handle Read / Write of H5 Flame files.

Created November 2016 by COOP team
"""

import os
import sys
import time
import logging
import subprocess
import commands
import copy
import numpy as np

from os.path import join, isdir
from math import pi
from h5py import File
from textwrap import dedent
from StringIO import StringIO
from tempfile import TemporaryFile

from geometry import (Parallelepiped, Sphere, Cylinder, Parallelogram, Disc,
                      ScatterShape2D, ScatterShape3D, NormalVector, Vector,
                      Point,)
from flamemetas import FlameMetas, VERSION

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
        mesh = File(avbp_mesh, 'r')
        sol = File(avbp_sol, 'r')
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
        if 'z' in mesh["Coordinates"].keys():
            z = mesh["Coordinates/z"].value
            z_in = z[above]
            self.log.debug("Generating 3D flame from AVBP. Bounding box :")
            self.log.debug(" ".join("{}".format(x) for x in [
                x_in.min(), x_in.max(),
                y_in.min(), y_in.max(),
                z_in.min(), z_in.max()]))
            self.shape = ScatterShape3D(
                    x_in.min(), x_in.max(),
                    y_in.min(), y_in.max(),
                    z_in.min(), z_in.max())
            self.metas.generation_method = "avbp_scalar_threshold_3D"
        else:
            self.shape = ScatterShape2D(
                    x_in.min(), x_in.max(),
                    y_in.min(), y_in.max())
            self.log.debug("Generating 2D flame from AVBP. Bounding box :")
            self.log.debug(" ".join("{}".format(x) for x in [
                x_in.min(), x_in.max(),
                y_in.min(), y_in.max()]))
            self.metas.generation_method = "avbp_scalar_threshold_2D"
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
        self.shape.inside_points = File(
                "dummy_flame.sol.h5", 'r')['Additionals/n_tau_flag'].value
        if not self.debug:
            os.remove("dummy_avbp.h5")
            os.remove("dummy_flametrans.h5")
            os.remove("dummy_flame.sol.h5")
            os.remove("dummy_flame.mesh.h5")
            os.remove("dummy_flame.asciiBound")

    def define_flame_disc(self, center, radius):
        """Define flame as parallelogram"""
        self.metas.generation_method = "analytical2D_disc"
        self.shape = Disc(center, radius)
        #self.metas.shape_params = np.hstack((center, radius))
        self.make_mesh()

    def define_flame_parallelogram(self, xref, vec1, vec2):
        """Define flame as parallelogram"""
        self.metas.generation_method = "analytical2D_parallelogram"
        self.shape = Parallelogram(xref, vec1, vec2)
        #self.metas.shape_params = np.hstack((xref, vec1, vec2))
        self.make_mesh()

    def define_flame_sphere(self, center, radius):
        """Define spherical flame"""
        self.metas.generation_method = "analytical3D_sphere"
        self.shape = Sphere(center, radius)
        #self.metas.shape_params = np.hstack((center, radius))
        self.make_mesh()

    def define_flame_cylinder(self, center, radius, vector):
        """Define cylindrical flame"""
        self.log.debug("Defining cylinder flame")
        self.metas.generation_method = "analytical3D_cylinder"
        self.shape = Cylinder(center, radius, vector)
        #self.metas.shape_params = np.hstack((center, radius, vector))
        self.make_mesh()

    def define_flame_parallelepiped(self, xref, vec1, vec2, vec3):
        """Define flame as parallelepiped"""
        self.metas.generation_method = "analytical3D_parallelepiped"
        self.shape = Parallelepiped(xref, vec1, vec2, vec3)
        #self.metas.shape_params = np.hstack((xref, vec1, vec2, vec3))
        self.make_mesh()

    def make_mesh(self):
        """Use bounding box to create box mesh containing flame geometry"""
        self.metas.ndim = self.shape.ndim
        self.metas.pt_min = copy.deepcopy(self.shape.vects["pt_min"])
        self.metas.pt_max = copy.deepcopy(self.shape.vects["pt_max"])
        self.exec_hip(self._get_hip_script_generate())
        self.read_meshpoints()

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

    def transform(self, trans, args):
        """Apply transformation to flame"""
        if trans == "t":
            self.shape.translate(args)
            self.metas.transform("translate", args)
            self.make_mesh()
        elif trans == "s":
            self.shape.scale(args)
            self.metas.transform("scale", args)
            self.make_mesh()
        elif trans == "r":
            # In this special case, flame generation method is lost (becomes ScatterShape)
            self.metas.transform("rotate", args)
            self._rotate(*args)
        else:
            raise ValueError

    def _rotate(self, axis, angle, degrees=True):
        """Write and rotate flame in hip, then recreate mesh"""
        self.write_h5(with_mesh=True)
        if not degrees: angle *= 180./np.pi
        rot = "tr ro {0} {1}".format(axis, angle)
        rotate_script = dedent("""\
          re -a flame_{0}.mesh.h5 -s flame_{0}.sol.h5
          {1}
          wr hd ./flame_{0}""".format(self.metas.name, rot))
        self.exec_hip(rotate_script)

        self.read_shape()

        interp_script = ["re -a flame_{0}.mesh.h5 -s -s flame_{0}.sol.h5"]
        interp_script += self._get_hip_script_generate().split("\n")
        interp_script.insert(3, "in gr 1")
        self.exec_hip('\n'.join(interp_script))
        self.make_mesh()

    def read_shape(self):
        """Read shape from meshpoints (only creates ScatterShape type)"""
        self.read_meshpoints()
        pt_min = self.meshpoints.min(axis=0)
        pt_max = self.meshpoints.max(axis=0)
        self.shape = (
                ScatterShape2D(pt_min[0], pt_max[0],
                               pt_min[1], pt_max[1]) if len(pt_min) == 2 else
                ScatterShape3D(pt_min[0], pt_max[0],
                               pt_min[1], pt_max[1],
                               pt_min[2], pt_max[2]))
            
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

    def read_meshpoints(self):
        """Read hdf5 meshfile from hip and store coordinates"""
        assert self.metas.ndim is not None, "ndim must be defined"
        with File(self.mesh_file, 'r') as f:
            if self.metas.ndim == 2:
                self.meshpoints = np.array((f['/Coordinates/x'].value,
                                            f['/Coordinates/y'].value)).T
            else:
                self.meshpoints = np.array((f['/Coordinates/x'].value,
                                            f['/Coordinates/y'].value,
                                            f['/Coordinates/z'].value)).T
            self.log.debug("Mesh data stored with {0} points".format(self.meshpoints[:,0].size))

    def write_h5(self, with_metas=True, with_mesh=False):
        """Write h5 file for flame"""
        assert self.meshpoints is not None
        if self.shape is None: self.read_shape()
        with File(self.flame_file, 'w') as f:
            f.create_group("/Additionals")
            if with_metas: self.metas.write_h5(f["/Additionals"].attrs)
            if 'analytical' in self.metas.generation_method:
                f['/Additionals/n_tau_flag'] = 1.0*self.shape.is_inside(self.meshpoints)
            else:
                f['/Additionals/n_tau_flag'] = self.shape.inside_points
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
        self.read_shape(path)
        self.make_mesh()
        self.read_meshpoints()

    def load_metas(self, path):
        """Load only the metas of flame at path"""
        self.log.debug("Loading metas for file " + path)
        with File(path, 'r') as flame:
            self.metas.load(flame['Additionals'].attrs)
        if float(self.metas.version) > float(VERSION):
            self.log.error("Cannot read future version flame. Please update your FlameTransfer app")
            raise ValueError

    #def recreate_shape(self, path):
    #    """Create self.shape again based on metas"""
    #    self.log.info("Recreating shape " + self.metas.generation_method)
    #    if self.metas.generation_method == "analytical2D_disc":
    #        center = self.metas.shape_params[:2]
    #        radius = self.metas.shape_params[2]
    #        self.shape = Disc(center, radius)
    #    elif self.metas.generation_method == "analytical2D_parallelogram":
    #        xref = self.metas.shape_params[:2]
    #        vec1 = self.metas.shape_params[2:4]
    #        vec2 = self.metas.shape_params[4:]
    #        self.shape = Parallelogram(xref, vec1, vec2)
    #    elif self.metas.generation_method == "analytical3D_sphere":
    #        center = self.metas.shape_params[:3]
    #        radius = self.metas.shape_params[3]
    #        self.shape = Sphere(center, radius)
    #    elif self.metas.generation_method == "analytical3D_cylinder":
    #        center = self.metas.shape_params[:3]
    #        radius = self.metas.shape_params[3]
    #        vector = self.metas.shape_params[4:]
    #        self.shape = Cylinder(center, radius, vector)
    #    elif self.metas.generation_method == "analytical3D_parallelepiped":
    #        xref = self.metas.shape_params[:3]
    #        vec1 = self.metas.shape_params[3:6]
    #        vec2 = self.metas.shape_params[6:9]
    #        vec3 = self.metas.shape_params[9:]
    #        self.shape = Parallelepiped(xref, vec1, vec2, vec3)
    #    elif "avbp_scalar_threshold" in self.metas.generation_method:
    #        if self.metas.ndim == 2:
    #            self.shape = ScatterShape2D(
    #                    self.metas.xmin, self.metas.xmax,
    #                    self.metas.ymin, self.metas.ymax,
    #                    )
    #        else:
    #            self.shape = ScatterShape3D(
    #                    self.metas.xmin, self.metas.xmax,
    #                    self.metas.ymin, self.metas.ymax,
    #                    self.metas.zmin, self.metas.zmax,
    #                    )
    #        with File(path, 'r') as f:
    #            self.shape.inside_points = f['Additionals/n_tau_flag'].value
    #    else:
    #        self.log.error("Unknown flame generation method")
    #        raise ValueError


