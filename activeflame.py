#!/usr/bin/env python
"""
ActiveFlame class to handle all information pertaining to an active flame.
Handle Read / Write of H5 Flame files.

Created November 2016 by COOP team
"""

import logging
import subprocess
import numpy as np

from os.path import join, isdir
from math import pi
from h5py import File
from textwrap import dedent

from geometry import Parallelepiped, Sphere, Cylinder, Parallelogram, Circle

GRID_SIZE = 32
VERSION = "0.1"

class ActiveFlame(object):
    """Flame holder class associated to a flame hdf5 file"""
    def __init__(self, name, **kwargs):
        self.log = logging.getLogger(__name__)
        self.debug = True
        self.meshpoints = None
        self.flame_file = None
        self.shape = None
        self.metas = FlameMetas(name=name)
        self.do_hip = False
        self.__dict__.update(kwargs)

    @staticmethod
    def compute_n_crocco(n1, area, p_mean, gamma):
        """Define gain (N2) using Crocco's analytical formulation"""
        return gamma * p_mean * area / (gamma - 1.) * n1

    def set_n1_tau(self, freq, tau, n1, area, p_mean, gamma):
        """Set N and tau from Crocco model (N1)"""
        self.set_n2_tau(freq, tau, compute_n_crocco(n1, area, p_mean, gamma))

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

    def define_hrthreshold_flame(self, avbp_mesh, avbp_sol, threshold):
        """Define flame from avbp Heat Release threshold"""
        raise NotImplementedError
        # TODO
        # Identify bounding box in source (hard!)
        # self.make_mesh()
        self.metas.generation_method = "scalar_and_threshold"
        # interpolate threshold on mesh
        # self.metas.volume = ??

    def define_flame_circle(self, center, radius):
        """Define flame as parallelogram"""
        self.metas.ndim = 2
        self.metas.generation_method = "analytical2D_circle"
        self.shape = Circle(center, radius)
        self.metas.shape_params = np.hstack((center, radius))
        self.update_metas(**self.shape.bounding_box())
        self.make_mesh_script()

    def define_flame_parallelogram(self, xref, vec1, vec2):
        """Define flame as parallelogram"""
        self.metas.ndim = 2
        self.metas.generation_method = "analytical2D_parallelogram"
        self.shape = Parallelogram(xref, vec1, vec2)
        self.metas.shape_params = np.hstack((xref, vec1, vec2))
        self.update_metas(**self.shape.bounding_box())
        self.make_mesh_script()

    def define_flame_sphere(self, center, radius):
        """Define spherical flame"""
        self.metas.ndim = 3
        self.metas.generation_method = "analytical3D_sphere"
        self.shape = Sphere(center, radius)
        self.metas.shape_params = np.hstack((center, radius))
        self.update_metas(**self.shape.bounding_box())
        self.make_mesh_script()

    def define_flame_cylinder(self, center, radius, vector):
        """Define cylindrical flame"""
        self.log.debug("Defining cylinder flame")
        self.metas.ndim = 3
        self.metas.generation_method = "analytical3D_cylinder"
        self.shape = Cylinder(center, radius, vector)
        self.metas.shape_params = np.hstack((center, radius, vector))
        self.update_metas(**self.shape.bounding_box())
        self.make_mesh_script()

    def make_mesh_script(self):
        """Use self.shape to create box mesh containing flame geometry"""
        with open('script.hip', 'w') as f:
            f.write(self._get_hip_script_generate())
        self.do_hip = True
        self.mesh_file = 'mesh_{0}.mesh.h5'.format(self.metas.name)

    def read_meshpoints(self):
        with File(self.mesh_file, 'r') as f:
            if self.metas.ndim == 2:
                self.meshpoints = np.array((f['/Coordinates/x'].value,
                                            f['/Coordinates/y'].value)).T
            elif self.metas.ndim == 3:
                self.meshpoints = np.array((f['/Coordinates/x'].value,
                                            f['/Coordinates/y'].value,
                                            f['/Coordinates/z'].value)).T
            self.log.debug("Mesh written, with {0} points".format(self.meshpoints[:,0].size))

    def _get_hip_script_generate(self):
        """Generate hip script for mesh generation"""
        assert self.metas.xmin is not None
        if self.metas.ndim == 2:
            hip_script = dedent("""\
              ge {0.xmin} {0.ymin} {0.xmax} {0.ymax} {0.grid_size} {0.grid_size}
              wr hd ./mesh_{0.name}
              bye""").format(self.metas)
        elif self.metas.ndim == 3:
            hip_script = dedent("""\
              ge {0.xmin} {0.ymin} {0.xmax} {0.ymax} {0.grid_size} {0.grid_size}
              co 3d {0.zmin} {0.zmax} {1} z
              wr hd ./mesh_{0.name}
              qu""").format(self.metas, self.metas.grid_size - 1)
        return hip_script

    def _get_hip_script_interpolate(self, src_mesh, src_sol, dest_mesh):
        """Generate hip script for flame interpolation on new mesh"""
        hip_script = dedent("""\
          re hd {src_mesh} {src_sol}
          re hd {dest_mesh}
          in gr 1
          wr hd -s ./sol_{self.metas.name}
          bye""").format(**locals())
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
            f['/Additionals/n_tau_flag'] = 1.0*self.shape.is_inside(self.meshpoints)
            f.create_group("/Parameters")
            f['/Parameters/ndim'] = np.array([self.metas.ndim])
            f['/Parameters/nnode'] = np.array([self.meshpoints.size])
            f['/Parameters/versionstring'] = "C3Sm_Flame"
        return path

    def write_n_tau(self, path):
        """Write ascii n_tau file for flame"""
        assert self.metas.n2_tau is not None
        path = 'n_tau_{}.dat'.format(self.metas.name)
        with open(path, 'w') as f:
            f.write("# Frequency N2 tau  - Generated using C3Sm\n")
            np.savetxt(f, self.metas.n2_tau.T)

    def load(self, path):
        """Load a flame at path"""
        self.load_metas(path)

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

    def recreate_shape(self):
        """Create self.shape again based on metas"""
        if self.metas.generation_method == "analytical2D_circle":
            center = self.metas.shape_params[:2]
            radius = self.metas.shape_params[2]
            self.shape = Circle(center, radius)
        if self.metas.generation_method == "analytical2D_parallelogram":
            xref = self.metas.shape_params[:2]
            vec1 = self.metas.shape_params[2:4]
            vec2 = self.metas.shape_params[4:]
            self.shape = Parallelogram(xref, vec1, vec2)
        if self.metas.generation_method == "analytical3D_sphere":
            center = self.metas.shape_params[:3]
            radius = self.metas.shape_params[3]
            self.shape = Sphere(center, radius)
        if self.metas.generation_method == "analytical3D_cylinder":
            center = self.metas.shape_params[:3]
            radius = self.metas.shape_params[3]
            vector = self.metas.shape_params[4:]
            self.shape = Cylinder(center, radius, vector)


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
