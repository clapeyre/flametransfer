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
import numpy as np

from os.path import join, isdir
from math import pi
from h5py import File
from textwrap import dedent
from StringIO import StringIO
from tempfile import TemporaryFile

from geometry import Parallelepiped, Sphere, Cylinder, Parallelogram, Circle, ScatterShape2D, ScatterShape3D

GRID_SIZE = 32
VERSION = "0.1"

class ActiveFlame(object):
    """Flame holder class associated to a flame hdf5 file"""
    def __init__(self, name, hip_exec, **kwargs):
        self.log = logging.getLogger(__name__)
        self.debug = True
        self.meshpoints = None
        self.flame_file = None
        self.shape = None
        self.metas = FlameMetas(name=name)
        self.hip_exec = hip_exec
        self.__dict__.update(kwargs)

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
                    sys.stdout.write(lines)
                    sys.stdout.flush()
            sys.stdout.write(output.read())
            sys.stdout.flush()
        print

    def compute_n_crocco(self, n1, area, p_mean, gamma):
        """Define gain (N2) using Crocco's analytical formulation"""
        return gamma * p_mean * area / (gamma - 1.) * n1

    def set_n1_tau(self, freq, tau, n1, area, p_mean, gamma):
        """Set N and tau from Crocco model (N1)"""
        self.set_n2_tau(freq, tau, self.compute_n_crocco(n1, area, p_mean, gamma))

    def set_n2_tau(self, freq, tau, n2):
        """Set N and tau from global values (same as AVSP internal)"""
        self.metas.n2_tau = np.array((freq, n2, tau), ndmin=2)
        if self.metas.n2_tau.shape[0] != 3:    # Ugly!! There's probably a better solution to keep orientation consistent
            self.metas.n2_tau = self.metas.n2_tau.T

    def set_n3_tau(self, freq, tau, n3, u_bar, q_bar):
        """Set N and tau from non dimensional values (N3)"""
        self.update_metas(**{"u_bar": u_bar, "q_bar": q_bar})
        self.set_n2_tau(freq, tau, n3 * u_bar / q_bar)

    def get_n3(self):
        """Get value of N3 using u_bar and q_bar"""
        if self.metas.q_bar is not None:
            return self.n2 * self.metas.q_bar / self.metas.u_bar
        else:
            self.log.error("q_bar info not available, cannot compute n3")

    def define_threshold_flame(self, avbp_mesh, avbp_sol, field, thresh):
        """Define flame from avbp scalar and threshold"""
        self.metas.generation_method = "avbp_scalar_threshold"
        mesh = File(avbp_mesh, 'r')
        sol = File(avbp_sol, 'r')
        def find_field(name):
            if field in name: return name
        data = sol[sol.visit(find_field)].value
        above = data > np.float(thresh)
        x = mesh["Coordinates/x"].value
        y = mesh["Coordinates/y"].value
        x_in = x[above]
        y_in = y[above]
        if 'z' in mesh["Coordinates"].keys():
            z = mesh["Coordinates/z"].value
            z_in = z[above]
            print x_in.min(), x_in.max(), y_in.min(), y_in.max(), z_in.min(), z_in.max()
            self.shape = ScatterShape3D(
                    x_in.min(), x_in.max(),
                    y_in.min(), y_in.max(),
                    z_in.min(), z_in.max())
        else:
            self.shape = ScatterShape2D(
                    x_in.min(), x_in.max(),
                    y_in.min(), y_in.max())
        self.make_mesh()
        with File("dummy_avbp.h5", 'w') as f:
            f.create_group
            f.create_group("/Additionals")
            f['/Additionals/n_tau_flag'] = 1.0*above
            f.create_group("/Parameters")
            f['/Parameters/ndim'] = np.array([self.metas.ndim])
            f['/Parameters/nnode'] = np.array([self.meshpoints.size])
            f['/Parameters/versionstring'] = "C3Sm_Flame"
        script = [
                "re hd -a {0} -s dummy_avbp.h5".format(avbp_mesh),
                "re hd {0}".format(self.mesh_file),
                "in gr 1",
                "wr hd -s dummy_flame",
                "qu",
                ]
        self.exec_hip('\n'.join(script))
        self.shape.inside_points = File(
                "dummy_flame.sol.h5", 'r')['Additionals/n_tau_flag'].value

    def define_flame_circle(self, center, radius):
        """Define flame as parallelogram"""
        self.metas.generation_method = "analytical2D_circle"
        self.shape = Circle(center, radius)
        self.metas.shape_params = np.hstack((center, radius))
        self.make_mesh()

    def define_flame_parallelogram(self, xref, vec1, vec2):
        """Define flame as parallelogram"""
        self.metas.generation_method = "analytical2D_parallelogram"
        self.shape = Parallelogram(xref, vec1, vec2)
        self.metas.shape_params = np.hstack((xref, vec1, vec2))
        self.make_mesh()

    def define_flame_sphere(self, center, radius):
        """Define spherical flame"""
        self.metas.generation_method = "analytical3D_sphere"
        self.shape = Sphere(center, radius)
        self.metas.shape_params = np.hstack((center, radius))
        self.make_mesh()

    def define_flame_cylinder(self, center, radius, vector):
        """Define cylindrical flame"""
        self.log.debug("Defining cylinder flame")
        self.metas.generation_method = "analytical3D_cylinder"
        self.shape = Cylinder(center, radius, vector)
        self.metas.shape_params = np.hstack((center, radius, vector))
        self.make_mesh()

    def define_flame_parallelepiped(self, xref, vec1, vec2, vec3):
        """Define flame as parallelepiped"""
        self.metas.generation_method = "analytical3D_parallelepiped"
        self.shape = Parallelogram(xref, vec1, vec2, vec3)
        self.metas.shape_params = np.hstack((xref, vec1, vec2, vec3))
        self.make_mesh()

    def make_mesh(self):
        """Use bounding box to create box mesh containing flame geometry"""
        self.update_metas(**self.shape.bounding_box())
        self.metas.ndim = self.shape.ndim
        self.mesh_file = 'mesh_{0}.mesh.h5'.format(self.metas.name)
        self.exec_hip(self._get_hip_script_generate())
        self.read_meshpoints()

    def _get_hip_script_generate(self):
        """Generate hip script for mesh generation"""
        assert self.metas.xmin is not None
        if self.metas.ndim == 2:
            hip_script = dedent("""\
              ge {0.xmin} {0.ymin} {0.xmax} {0.ymax} {0.grid_size} {0.grid_size}
              wr hd ./mesh_{0.name}
              qu
              """).format(self.metas)
        elif self.metas.ndim == 3:
            hip_script = dedent("""\
              ge {0.xmin} {0.ymin} {0.xmax} {0.ymax} {0.grid_size} {0.grid_size}
              co 3d {0.zmin} {0.zmax} {1} z
              wr hd ./mesh_{0.name}
              qu
              """).format(self.metas, self.metas.grid_size - 1)
        return hip_script

    def read_meshpoints(self):
        """Read hdf5 meshfile from hip and store coordinates"""
        with File(self.mesh_file, 'r') as f:
            if self.metas.ndim == 2:
                self.meshpoints = np.array((f['/Coordinates/x'].value,
                                            f['/Coordinates/y'].value)).T
            elif self.metas.ndim == 3:
                self.meshpoints = np.array((f['/Coordinates/x'].value,
                                            f['/Coordinates/y'].value,
                                            f['/Coordinates/z'].value)).T
            self.log.debug("Mesh data stored with {0} points".format(self.meshpoints[:,0].size))

    def _get_hip_script_interpolate(self, src_mesh, src_sol, dest_mesh):
        """Generate hip script for flame interpolation on new mesh"""
        hip_script = dedent("""\
          re hd {src_mesh} {src_sol}
          re hd {dest_mesh}
          in gr 1
          wr hd -s ./sol_{self.metas.name}
          qu""").format(**locals())
        return hip_script

    def write_h5(self, with_metas=True):
        """Write h5 file for flame"""
        assert self.meshpoints is not None
        if self.shape is None: self.recreate_shape()
        path = "flame_{0}.h5".format(self.metas.name)
        self.flame_file = path
        with File(path, 'w') as f:
            f.create_group("/Additionals")
            if with_metas: self.metas.write(f["/Additionals"].attrs)
            if 'analytical' in self.metas.generation_method:
                f['/Additionals/n_tau_flag'] = 1.0*self.shape.is_inside(self.meshpoints)
            else:
                f['/Additionals/n_tau_flag'] = self.shape.inside_points
            f.create_group("/Parameters")
            f['/Parameters/ndim'] = np.array([self.metas.ndim])
            f['/Parameters/nnode'] = np.array([self.meshpoints.size])
            f['/Parameters/versionstring'] = "C3Sm_Flame"
        return path

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
        self.recreate_shape(path)
        self.make_mesh()
        self.read_meshpoints()

    def load_metas(self, path):
        """Load only the metas of flame at path"""
        flame = File(path, 'r')
        self.metas.load(flame['Additionals'].attrs)
        if float(self.metas.version) > float(VERSION):
            self.log.error("Cannot read future version flame. Please update your FlameTransfer app")
            raise ValueError

    def update_metas(self, **kwargs):
        """Update meta data"""
        self.metas.__dict__.update(**kwargs)

    def recreate_shape(self, path):
        """Create self.shape again based on metas"""
        if self.metas.generation_method == "analytical2D_circle":
            center = self.metas.shape_params[:2]
            radius = self.metas.shape_params[2]
            self.shape = Circle(center, radius)
        elif self.metas.generation_method == "analytical2D_parallelogram":
            xref = self.metas.shape_params[:2]
            vec1 = self.metas.shape_params[2:4]
            vec2 = self.metas.shape_params[4:]
            self.shape = Parallelogram(xref, vec1, vec2)
        elif self.metas.generation_method == "analytical3D_sphere":
            center = self.metas.shape_params[:3]
            radius = self.metas.shape_params[3]
            self.shape = Sphere(center, radius)
        elif self.metas.generation_method == "analytical3D_cylinder":
            center = self.metas.shape_params[:3]
            radius = self.metas.shape_params[3]
            vector = self.metas.shape_params[4:]
            self.shape = Cylinder(center, radius, vector)
        elif self.metas.generation_method == "analytical3D_parallelepiped":
            xref = self.metas.shape_params[:3]
            vec1 = self.metas.shape_params[3:6]
            vec2 = self.metas.shape_params[6:9]
            vec3 = self.metas.shape_params[9:]
            self.shape = Parallelepiped(xref, vec1, vec2, vec3)
        elif self.metas.generation_method == "avbp_scalar_threshold":
            if self.metas.ndim == 2:
                self.shape = ScatterShape2D(
                        self.metas.xmin, self.metas.xmax,
                        self.metas.ymin, self.metas.ymax,
                        )
            else:
                self.shape = ScatterShape3D(
                        self.metas.xmin, self.metas.xmax,
                        self.metas.ymin, self.metas.ymax,
                        self.metas.zmin, self.metas.zmax,
                        )
            self.shape.inside_points = File(path, 'r')['Additionals/n_tau_flag'].value


class FlameMetas(object):
    """Meta data associated to a flame"""
    def __init__(self, **kwargs):
        self.name = None
        self.ndim = None
        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.zmin = None
        self.zmax = None
        self.ref_point = None
        self.ref_vect = None
        self.geom_ref_point = None
        self.geom_ref_vect = None
        self.n2_tau = None
        self.q_bar = None
        self.u_bar = None
        self.shape_params = None
        self.generation_method = None
        self.grid_size = GRID_SIZE
        self.version = VERSION
        self.__dict__.update(**kwargs)

    def write(self, attribute):
        """Write to hdf5 attributes and check mandatories are filled"""
        mandatory_attrs = ("xmin xmax ymin ymax ref_point ref_vect n2_tau" 
                           "generation_method").split()
        for key, value in self.__dict__.iteritems():
            if key in mandatory_attrs: assert value is not None
            if value is not None: attribute[key] = value

    def load(self, attribute):
        """Load self from hdf5 attribute"""
        for key, value in attribute.iteritems():
            setattr(self, key, value)
