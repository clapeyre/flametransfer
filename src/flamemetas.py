#!/usr/bin/env python
"""
FlameMetas class to handle all meta information pertaining to an active flame.

Created January 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""

import json

from collections import OrderedDict
from numpy import array as nparray

import geometry
from constants import VERSION

class FlameMetas(object):
    """Meta data associated to a flame"""
    def __init__(self, **kwargs):
        # Override of __setattr__ prevents declarations here. 
        # Calling object.__setattr__ is an approved workaround 
        object.__setattr__(self, "mandatory_vals",
                           "name version grid_size ndim generation_method "
                           "n2_tau".split())
        object.__setattr__(self, "mandatory_vecs", 
                           "pt_min pt_max pt_ref vec_ref".split())
        object.__setattr__(self, "static", OrderedDict((
                             (val, None) for val in self.mandatory_vals)))
        object.__setattr__(self, "vects", OrderedDict((
                             (vec, None) for vec in self.mandatory_vecs)))
        self.static["version"] = VERSION
        self.static["grid_size"] = GRID_SIZE
        self.static.update(**kwargs)

    def __setattr__(self, name, value):
        """Override set attribute behavior
        
        Everything goes to static, except if already in self.vects
        See set_vect to ensure value goes to self.vects
        """
        if name in self.vects.keys():
            self.vects[name] = value
        else:
            self.static[name] = value

    def __getattr__(self, name):
        """Override get attribute behavior
        
        Look for attributes in static and vects. None other are accepted
        """
        if name in self.static.keys():
            return self.static[name]
        elif name in self.vects.keys():
            if self.vects[name] is None:
                return self.vects[name]
            else:
                return self.vects[name]()
        else:
            raise AttributeError(name)

    def set_vect(self, name, vect):
        """Set new attribute ensured to be a vector"""
        self.vects[name] = vect

    def write_h5(self, attribute):
        """Check mandatories are filled and write to hdf5 attribute"""
        for key, value in self.get_metas().iteritems():
            if key in self.mandatory_vals + self.mandatory_vecs:
                assert value is not None, "{} is mandatory for write".format(key)
            attribute[key] = value

    def load(self, attribute):
        """Load self from hdf5 attribute"""
        self.static.update(((key, attribute[key]) for key in attribute["static"]))
        self.vects.update((
                (key, getattr(geometry, self.vect_types[i])(attribute[k]))
                for i,k in enumerate(attribute["vects"])))

    def get_metas(self):
        """Get all the meta data for the flame"""
        out = self.static.copy()
        out.update({
            vect: self.vects[vect]() if self.vects[vect] is not None else None
            for vect in self.vects.keys()})
        out["static"] = self.static.keys()
        out["vects"] = self.vects.keys()
        out["vect_types"] = [v.__class__.__name__ for k,v in self.vects.items()]
        return out

    def write_json(self, filehandle):
        """Dump all static as JSON file
        
        Twist: all numpy arrays must be converted to lists before JSON dump
        The `numpies` field is used to know who was a numpy before dump (for
        reconstruction).
        """
        import json
        out = self.get_metas()
        out['numpies'] = out["vects"][:]
        for k, v in out.iteritems():
            try:
                out[k] = v.tolist()
                out['numpies'].append(k)
            except AttributeError:
                pass
        json.dump(out, filehandle, indent=4, separators=(',', ': '))

    def transform(self, method, *args):
        """Apply transformation to metas"""
        [getattr(v, method)(*args) if v is not None else None
         for v in self.vects.values()]
