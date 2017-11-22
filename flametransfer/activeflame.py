#!/usr/bin/env python
"""
ActiveFlame class to handle all information pertaining to an active flame.
Handle Read / Write of H5 Flame files.

Created November 2016 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""
import os
import logging
import copy
import shutil

from os.path import isfile, realpath, dirname, join
from textwrap import dedent

import numpy as np

from h5py import File

from geometry import (Rectangle, Sphere, Cylinder, Brick, Disc, ScatterShape2D,
                      ScatterShape3D, Point, json2shape, shape2json)
from flamemetas import FlameMetas
from constants import version_checker, DEBUG
from tools import visu


class ActiveFlame(object):
    """Flame holder class associated to a flame hdf5 file"""
    def __init__(self, name, hip_wrapper, **kwargs):
        self.log = logging.getLogger(__name__)
        self.meshpoints = None
        self.inside_points = None
        self.shape = None
        self.metas = FlameMetas(name=name)
        self.hip_wrapper = hip_wrapper
        self.__dict__.update(kwargs)

    def __deepcopy__(self, memo):
        newone = type(self)(self.name, self.hip_wrapper)
        newone.__dict__.update(self.__dict__)
        newone.metas = copy.deepcopy(self.metas, memo)
        return newone

    @property
    def name(self):
        """Flame name, stored in metas"""
        return self.metas.name

    @property
    def mesh_file(self):
        """Mesh file name"""
        return '{0}.mesh.h5'.format(self.name)

    @property
    def mesh_file_hip(self):
        """Mesh file name as written by hip"""
        return '{0}.mesh.h5'.format(self.name)

    @property
    def flame_file(self):
        """Flame file name"""
        return '{0}.flame.h5'.format(self.name)

    @property
    def flame_file_hip(self):
        """Flame file name as written by hip"""
        return '{0}.sol.h5'.format(self.name)

    @property
    def template(self):
        """Name of template file"""
        return join(dirname(realpath(__file__)), "..", "template",
                    "template{}d.mesh.h5".format(self.metas.ndim))

    def harmonize_flame_name(self):
        """Hip fixes the name of output file, not as we like"""
        if isfile(self.mesh_file_hip):
            pass  # Right now, hip writes the correct file name
        if isfile(self.flame_file_hip):
            os.rename(self.flame_file_hip, self.flame_file)
        # if isfile(self.mesh_file_hip.replace('h5', 'xmf')):
        #     with open(self.mesh_file_hip.replace('h5', 'xmf')) as xmf:
        #         xmf = xmf.readlines()
        #     xmf = [line.replace(self.flame_file_hip, self.flame_file)
        #            for line in xmf]
        #     with open(self.mesh_file_hip.replace('h5', 'xmf'), 'w') as out:
        #         out.writelines(xmf)

    def __getstate__(self):
        dic = dict(self.__dict__)
        del dic['log']
        return dic

    def __setstate__(self, dic):
        self.__dict__.update(dic)
        self.log = logging.getLogger(__name__)

    @staticmethod
    def compute_n_crocco(n_1, area, p_mean, gamma):
        """Define gain (N2) using Crocco's analytical formulation"""
        return gamma * p_mean * area / (gamma - 1.) * n_1

    def set_n1_tau(self, freq, n_1, tau, area, p_mean, gamma):
        """Set N and tau from Crocco model (N1)"""
        self.set_n2_tau(freq,
                        self.compute_n_crocco(n_1, area, p_mean, gamma),
                        tau)

    def set_n2_tau(self, freq, n_2, tau):
        """Set N and tau from global values (same as AVSP internal)"""
        self.metas.n2_tau = np.array((freq, n_2, tau), ndmin=2)
        # Ugly!! There's probably a better way to keep orientation consistent
        if self.metas.n2_tau.shape[0] != 3:
            self.metas.n2_tau = self.metas.n2_tau.T

    def set_n3_tau(self, freq, n_3, tau, u_bar, q_bar):
        """Set N and tau from non dimensional values (N3)"""
        self.metas.u_bar = u_bar
        self.metas.q_bar = q_bar
        self.set_n2_tau(freq, n_3 * u_bar / q_bar, tau)

    def get_n3(self):
        """Get value of N3 using u_bar and q_bar"""
        if hasattr(self.metas, "q_bar"):
            return (self.metas.n2_tau[1, :] * self.metas.q_bar /
                    self.metas.u_bar)
        else:
            self.log.error("q_bar info not available, cannot compute n3")

    def define_flame_avbp_threshold(self, avbp_mesh, avbp_sol, field, thresh):
        """Define flame from avbp scalar and threshold"""
        self.log.debug("Generating scattershape from AVBP")
        with File(avbp_sol, 'r') as sol:
            def find_field(name):
                """h5py method for finding fields"""
                if field.lower() in name.lower():
                    return name
            selected_field = sol.visit(find_field)
            assert selected_field is not None, (
                "no field containing {0} found in {1}"
                .format(field, avbp_sol))
            data = sol[sol.visit(find_field)].value
            flame_flag = data > np.float(thresh)
        self.metas.avbp_mesh = avbp_mesh
        self.metas.avbp_sol = avbp_sol
        self.metas.field = selected_field
        self.metas.threshold = thresh
        self.define_flame_scatter(avbp_mesh, flame_flag)

    def define_flame_scatter(self, input_mesh, flame_flag):
        """Use flame_flag on input_mesh to define ScatterPoint flame"""
        with File(input_mesh, 'r') as mesh:
            x_flame = mesh["Coordinates/x"].value[flame_flag]
            y_flame = mesh["Coordinates/y"].value[flame_flag]
            pt_min = [x_flame.min(), y_flame.min()]
            pt_max = [x_flame.max(), y_flame.max()]
            klass = ScatterShape2D
            self.metas.generation_method = "avbp_scalar_threshold_2D"
            if 'z' in mesh["Coordinates"].keys():
                self.metas.generation_method = "avbp_scalar_threshold_3D"
                z_flame = mesh["Coordinates/z"].value[flame_flag]
                pt_min += [z_flame.min()]
                pt_max += [z_flame.max()]
                klass = ScatterShape3D
        self._shape_metas(klass(pt_min, pt_max))
        self.make_mesh()
        version = "flametransfer_v" + self.metas.version
        with File("dummy.h5", 'w') as flame:
            flame.create_group("/Average")
            flame['/Average/flames'] = 1.0*flame_flag
            flame.create_group("/Parameters")
            flame['/Parameters/ndim'] = np.array([self.metas.ndim])
            flame['/Parameters/nnode'] = np.array([self.meshpoints.shape[0]])
            flame['/Parameters/versionstring'] = version
        with File("dummy_flametrans.h5", 'w') as flame:
            flame.create_group("/Average")
            flame['/Average/flames'] = 0.0*flame_flag
            flame.create_group("/Parameters")
            flame['/Parameters/ndim'] = np.array([self.metas.ndim])
            flame['/Parameters/nnode'] = np.array([self.meshpoints.shape[0]])
            flame['/Parameters/versionstring'] = version
        script = [
            "re hd -a {0} -s dummy.h5".format(input_mesh),
            "re hd -a {0} -s dummy_flametrans.h5".format(self.mesh_file),
            "set in-rim 0.1",
            "in gr 1",
            "wr hd dummy_flame",
            "qu",
            ]
        self.hip_wrapper.execute('\n'.join(script))
        with File("dummy_flame.sol.h5", 'r') as flame:
            self.inside_points = flame["/Average/flames"].value
        if not DEBUG:
            os.remove("dummy.h5")
            os.remove("dummy_flametrans.h5")
            os.remove("dummy_flame.sol.h5")
            os.remove("dummy_flame.mesh.h5")
            os.remove("dummy_flame.asciiBound")

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
        """Set metas that define the shape"""
        self.metas.ndim = shape.ndim
        self.metas.pt_min = copy.deepcopy(shape.vects.pt_min)
        self.metas.pt_max = copy.deepcopy(shape.vects.pt_max)
        self.metas.volume = shape.volume
        self.metas.shape_params = shape2json(shape)

    def _make_analytic_shape(self, shape):
        self._shape_metas(shape)
        self.make_meshpoints()
        self.inside_points = 1.0*shape.is_inside(self.meshpoints)

    def make_mesh(self):
        """Create and read mesh for flame

        Use bounding box to create box mesh containing flame geometry
        """
        assert self.metas.pt_min is not None, (
            "Bounding box needed for mesh gen")
        self.make_meshpoints()
        shutil.copy(self.template, self.mesh_file)
        with File(self.mesh_file, "a") as mesh:
            mesh["Coordinates/x"][...] = self.meshpoints[:, 0]
            mesh["Coordinates/y"][...] = self.meshpoints[:, 1]
            if self.metas.ndim == 3:
                mesh["Coordinates/z"][...] = self.meshpoints[:, 2]
        self.log.info("Wrote mesh file to " + self.mesh_file)

    def make_meshpoints(self):
        """Store coordinates as hip writes them"""
        assert self.metas.ndim in [2, 3], "ndim must be defined"
        pt_min = self.metas.pt_min
        pt_max = self.metas.pt_max
        size = pt_max - pt_min
        self.meshpoints = np.zeros((self.metas.grid_size**self.metas.ndim,
                                    self.metas.ndim))
        with File(self.template, "r") as templ:
            self.meshpoints[:, 0] = pt_min[0] + templ["Coordinates/x"]*size[0]
            self.meshpoints[:, 1] = pt_min[1] + templ["Coordinates/y"]*size[1]
            if self.metas.ndim == 3:
                self.meshpoints[:, 2] = (pt_min[2]
                                         + templ["Coordinates/z"]*size[2])
        self.log.debug("Mesh data stored with {0} points"
                       .format(self.meshpoints[:, 0].size))

    def _get_hip_script_generate(self):
        """Write hip script for mesh generation"""
        line_3d = ("co 3d {0[2]} {1[2]} {2} z".format(self.metas.pt_min,
                                                      self.metas.pt_max,
                                                      self.metas.grid_size - 1)
                   if self.metas.ndim == 3 else "")
        hip_script = dedent("""\
          ge {0[0]} {0[1]} {1[0]} {1[1]} {2} {2}
          {3} 
          """).format(self.metas.pt_min, self.metas.pt_max,
                      self.metas.grid_size, line_3d)
        return hip_script

    def transform(self, trans, *args):
        """Apply transformation to flame

        A "radical" choice here is to discard any analytic shape information if
        a transformation is applied. While a smarter approach could apply the
        transformation to the analytic shape parameters, this need has not yet
        been expressed.
        """
        shape = json2shape(self.metas.shape_params)
        if trans == "t":
            self.metas.transform("translate", *args)
            shape.translate(*args)
        elif trans == "s":
            self.metas.transform("scale", *args)
            shape.scale(*args)
        elif trans == "r":
            # Special case: box must be rotated, then interpolated back on a
            # cartesian grid aligned with x,y[,z]
            self._rotate(*args)
            # transform metas, and overwrite metas for min/max points
            self.metas.transform("rotate", *args)
            self.metas.pt_min = Point(self.meshpoints.min(0))
            self.metas.pt_max = Point(self.meshpoints.max(0))
            # tranform shape to recover correct json
            shape.rotate(*args)
            shape.vects.pt_min = Point(self.meshpoints.min(0))
            shape.vects.pt_max = Point(self.meshpoints.max(0))
        else:
            raise ValueError("Unknown transformation, should be in [t, s, r]")
        self.make_meshpoints()
        self._shape_metas(shape)

    def _convert2scatter(self):
        """Convert this flame type to ScatterShape"""
        klass = ScatterShape3D if self.metas.ndim == 3 else ScatterShape2D
        self._shape_metas(klass(self.metas.pt_min, self.metas.pt_max))

    def set_inside_pts(self):
        """Set points inside flame according to flame file"""
        with File(self.flame_file, 'r') as flame:
            self.inside_points = flame["/Average/flames"].value

    def read_meshpoints(self):
        """Read coordinates of the current mesh.

        This mesh might not match the FlameTransfer format (e.g. cartesian
        aligned) so self.meshpoints is not updated.
        """
        meshpoints = np.zeros((self.metas.grid_size**self.metas.ndim,
                               self.metas.ndim))
        with File(self.mesh_file, 'r') as mesh:
            meshpoints[:, 0] = mesh["/Coordinates/x"].value
            meshpoints[:, 1] = mesh["/Coordinates/y"].value
            if self.metas.ndim == 3:
                meshpoints[:, 2] = mesh["/Coordinates/z"].value
        return meshpoints

    def _rotate(self, axis, angle, degrees=True):
        """Write and rotate flame in hip, then recreate x,y[,z] aligned mesh"""
        self.make_mesh()
        self.write_h5()

        # Perform rotation of the mesh with hip
        if not degrees:
            angle *= 180./np.pi
        rotate_script = dedent("""\
          re hd -a {0.mesh_file} -s {0.flame_file}
          tr ro {1} {2}
          wr hd ./{3}
          """.format(self, axis, angle, self.name))
        self.hip_wrapper.execute(rotate_script)
        self.harmonize_flame_name()

        # Extract all flame point coordinates to create new cartesian aligned
        # bounding box
        meshpoints = self.read_meshpoints()
        self.set_inside_pts()
        flamepoints = meshpoints[self.inside_points > 0.01]
        self.metas.pt_min = Point(flamepoints.min(0))
        self.metas.pt_max = Point(flamepoints.max(0))
        self.make_meshpoints()

        # Interpolate the rotated flame on the new cartesian aligned grid
        interp_script = ["re hd -a {0.mesh_file} -s {0.flame_file}"
                         .format(self)]
        interp_script += self._get_hip_script_generate().split("\n")
        interp_script += ["se ch 0",
                          "se in-rim 0",
                          # WAITING FOR HIP INTERPOLATION BUG CORRECTION
                          # "se bc-ty * n",
                          # "se in-recoType flag",
                          "in gr 1",
                          "wr hd ./{0}".format(self.name)]
        self.hip_wrapper.execute('\n'.join(interp_script))
        self.harmonize_flame_name()
        self.set_inside_pts()
        self.write_h5()

    def export_avsp(self, avsp_mesh, avsp_sol, number):
        """Export flame to AVSP solution"""
        self.make_mesh()
        name = "flame_{}.h5".format(number)
        self.write_h5(number=number, name=name)
        # with File("avsp.sol.h5", "r") as f: print f["Average/flames"].value
        self.hip_wrapper.execute(self._get_expavsp_hip_script(
            self.mesh_file, name, avsp_mesh, avsp_sol, number))
        if not DEBUG:
            os.remove(name)

    @staticmethod
    def _get_expavsp_hip_script(src_mesh, src_sol, avsp_mesh, avsp_sol,
                                number):
        """Generate hip script for flame interpolation on AVSP mesh"""
        return dedent("""\
          se ch 0
          se in-recoType flag
          se in-rim 0
          se bc-ty * n
          re hd -a {0} -s {1}
          re hd -a {2} -s {3}
          in gr 1
          wr hd -a avsp_tmp{4}
          """).format(src_mesh, src_sol, avsp_mesh, avsp_sol, number)

    def write_h5(self, number=1, name=None):
        """Write h5 file for flame"""
        assert self.meshpoints is not None, "Flame not read/generated yet"
        name = name if name is not None else self.flame_file
        with File(name, 'w') as flame:
            flame.create_group("/Average")
            flame['/Average/flames'] = np.where(self.inside_points > 0.1,
                                                float(number), 0.)
            self.write_group(flame, number=number)
            flame.create_group("/Parameters")
            flame['/Parameters/ndim'] = np.array([self.metas.ndim])
            flame['/Parameters/nnode'] = np.array([self.meshpoints.shape[0]])
            flame['/Parameters/versionstring'] = ("flametransfer_v"
                                                  + self.metas.version)
            self.log.info("Wrote flame file to " + self.flame_file)

    def write_full(self):
        """Write mesh, flame and xmf file"""
        self.make_mesh()
        self.write_h5()
        output = visu(self.mesh_file, self.flame_file)

    def write_group(self, hdf_obj, number=1):
        """Write full flame in /Flames/00X of hdf file"""
        if "Flames" not in hdf_obj.keys():
            hdf_obj.create_group("Flames")
        grp = hdf_obj.create_group("/Flames/{:03}".format(number))
        grp["n2_tau"] = self.metas.n2_tau
        grp["pt_ref"] = self.metas.pt_ref
        grp["vec_ref"] = self.metas.vec_ref
        grp["flame"] = self.inside_points
        self.metas.write_h5(grp.attrs)

    def load_group(self, hdf_obj, number=1):
        """Load all metas for flame in /Flames/00X of hdf file"""
        assert "Flames" in hdf_obj.keys(), (
            "HDF object doesn't contain Flame data")
        grp = hdf_obj["/Flames/{:03}".format(number)]
        version_checker(grp.attrs['version'])
        self.metas.load(grp.attrs)

    def read_inside_pts(self, hdf_obj, number=1):
        """Read inside points in /Flames/number of hdf file"""
        self.inside_points = hdf_obj['/Flames/{:03}/flame'
                                     .format(number)].value

    def load(self, path, number=1):
        """Load a flame at path"""
        self.log.info("Loading metas for file " + path)
        with File(path, 'r') as flame:
            self._read_metas(flame, number=number)
            self.read_inside_pts(flame, number=number)
        self.log.info("Loaded flame named " + self.name)
        self.log.debug("Contents :\n" + "\n"
                       .join("{0}: {1}".format(k, v)
                             for k, v in self.metas.__dict__.items()))
        self.make_meshpoints()

    def load_metas(self, path, number=1):
        """Wrap internal _read_metas function for reading only metas"""
        with File(path, 'r') as flame:
            self._read_metas(flame, number=number)

    def _read_metas(self, flame, number=1):
        """Load the metas contained in an hdf group"""
        grp = flame["/Flames/{:03}".format(number)]
        version_checker(grp.attrs['version'])
        self.metas.load(grp.attrs)
