#!/usr/bin/env python
"""
FlameMetas class to handle all meta information pertaining to an active flame.

Created January 2017 by COOP team
"""

import json

from numpy import array as nparray

from geometry import Vector

GRID_SIZE = 32
VERSION = "17.1"

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
        self.ignored = ["ignored"]
        self.transforms = []
        self.grid_size = GRID_SIZE
        self.version = VERSION
        self.__dict__.update(**kwargs)

    def write_h5(self, attribute):
        """Write to hdf5 attribute and check mandatories are filled"""
        mandatory_attrs = ("xmin xmax ymin ymax ref_point ref_vect n2_tau" 
                           "generation_method").split()
        for key, value in self.__dict__.iteritems():
            if key in mandatory_attrs: assert value is not None
            elif key not in self.ignored and value is not None:
                if key == "transforms": value = json.dumps(value)
                attribute[key] = value

    def load(self, attribute):
        """Load self from hdf5 attribute"""
        for key, value in attribute.iteritems():
            if key in self.ignored: continue
            else:
                if "transforms" in key: value = json.loads(value)
                setattr(self, key, value)

    def apply_trans(self):
        """Apply all transformations to all points and vectors"""
        pt_min = [self.xmin, self.ymin]
        pt_max = [self.xmax, self.ymax]
        if self.ndim == 3:
            pt_min += [self.zmin]
            pt_max += [self.zmax]
        points = [pt_min, pt_max]
        if self.ref_point is not None: points += [self.ref_point]
        if self.geom_ref_point is not None: points += [self.geom_ref_point]
        pt_update = [Vector(self.ndim, *arr[:]) for arr in points]
        vectors = []
        if self.ref_vect is not None: vectors.append(self.ref_vect)
        if self.geom_ref_vect is not None: vectors.append(self.geom_ref_vect)
        vec_update = [Vector(self.ndim, *arr[:]) for arr in vectors]

        for trans, args in self.transforms:
            if trans == "tr":
                [p.translate(args) for p in pt_update]
            elif trans == "sc":
                [p.scale(args) for p in pt_update]
            elif trans == "ro":
                [p.rotate(*args) for p in pt_update]
                [p.rotate(*args) for p in vec_update]
            elif trans == "mi":
                axis_dict = {'x':1, 'y':2, 'z':3}
                [p.mirror(args) for p in pt_update]
                [p.mirror(args) for p in vec_update]

        if self.ndim == 2:
            self.tr_xmin, self.tr_ymin = pt_update[0].vect
            self.tr_xmax, self.tr_ymax = pt_update[1].vect
        else:
            self.tr_xmin, self.tr_ymin, self.tr_zmin = pt_update[0].vect
            self.tr_xmax, self.tr_ymax, self.tr_zmax = pt_update[1].vect
        if self.ref_point is not None:
            self.tr_ref_point = pt_update[2].vect
        if self.geom_ref_point is not None:
            self.tr_geom_ref_point = pt_update[3].vect
        if self.ref_vect is not None:
            self.tr_ref_vect = vec_update[0].vect
        if self.geom_ref_vect is not None:
            self.tr_geom_ref_vect = vec_update[1].vect

    def write_json(self, filehandle):
        """Dump all metas as JSON file
        
        Twist: all numpy arrays must be converted to list before JSON dump
        A `numpies` field is added to know who was a numpy before dump (for
        reconstruction).
        """
        import json
        numpies = []
        out = {}
        for k, v in self.__dict__.iteritems():
            try:
                out[k] = v.tolist()
                numpies.append(k)
            except AttributeError:
                out[k] = v
        out['numpies'] = numpies 
        json.dump(out, filehandle, indent=4, separators=(',', ': '))
