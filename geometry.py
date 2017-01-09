#/usr/bin/env python
# -*- coding: utf-8 -*-
"""geometry.py

TODO

Created Nov 2016 by COOP team
"""
import numpy as np
import math
from abc import abstractmethod

def intersect(start1, end1, start2, end2):
    """Given 4 points, find intersection of 2 lines. Works in any dimension"""
    den = np.linalg.det(np.array((end1 - start1, end2 - start2)))
    num = np.linalg.det(np.array((start2 - start1, end2 - start2)))
    return start1 + num / den * (end1 - start1)


class Vector(object):
    """Base vector type"""
    def __init__(self, vect):
        self.ndim = len(vect)
        self.vect = np.array(vect).reshape(-1, 1)
        self.length = np.linalg.norm(self.vect)

    def __call__(self): return self.vect.reshape(-1)

    def rotate(self, axis, angle, degrees=True):
        if degrees: angle *= np.pi/180.
        if self.ndim == 2:
            rotate = np.mat([[math.cos(angle), -math.sin(angle)],
                             [math.sin(angle), math.cos(angle)]])
        else:
            if axis == 'x':
                rotate = np.mat([[1, 0, 0],
                                 [0, math.cos(angle), -math.sin(angle)],
                                 [0, math.sin(angle), math.cos(angle)]])
            elif axis == 'y':
                rotate = np.mat([[math.cos(angle), 0, math.sin(angle)],
                                 [0, 1, 0],
                                 [-math.sin(angle), 0, math.cos(angle)]])
            elif axis == 'z':
                rotate = np.mat([[math.cos(angle), -math.sin(angle), 0],
                                 [math.sin(angle), math.cos(angle), 0],
                                 [0, 0, 1]])
        self.vect = rotate * self.vect[:]

    def scale(self, vect):
        """Scale (negative accepted, 0 refused)"""
        vect = np.array(vect)
        assert all(abs(vect) > 1e-15), "Scaling must not have 0 value"
        self.vect *= vect.reshape(-1, 1)

    def translate(self, vect):
        """Vectors are not affected by translation"""
        pass


class NormalVector(Vector):
    """More restrictive vector type"""
    def __init__(self, vect):
        Vector.__init__(self, vect)
        assert self.length > 1e-15, "Vector cannot have 0 length "
        self.vect /= self.length

    def scale(self, vect):
        """Normalized vectors are not scalable"""
        pass


class Point(Vector):
    def translate(self, vect):
        """Points are affected by translation"""
        self.vect += np.array(vect).reshape(-1, 1)


class Shape(object):
    def __init__(self, ndim):
        self.ndim = ndim
        self.vects = {}

    def rotate(self, axis, angle, degrees=True):
        [p.rotate(axis, angle, degrees) for p in self.vects.values()]

    def scale(self, vect):
        [p.scale(vect) for p in self.vects.values()]

    def translate(self, vect):
        [p.translate(vect) for p in self.vects.values()]


class Shape2D(Shape):
    def __init__(self):
        Shape.__init__(self, 2)

    def check(self):
        assert self.area > 1.e-12


class Shape3D(Shape):
    def __init__(self):
        Shape.__init__(self, 3)

    def check(self):
        assert self.volume > 1.e-18


class ScatterShape2D(Shape2D):
    def __init__(self, pt_min, pt_max):
        Shape2D.__init__(self)
        self.vects["pt_min"] = Point(pt_min)
        self.vects["pt_max"] = Point(pt_max)
        self.inside_points = None
        self.check()

    @property
    def area(self):
        """Get shape area"""
        # TODO: The real area should be computed...
        pt_min = self.vects["pt_min"]()
        pt_max = self.vects["pt_max"]()
        return np.prod(pt_max - p_min)
    volume = area


class ScatterShape3D(Shape3D):
    def __init__(self, pt_min, pt_max):
        Shape3D.__init__(self)
        self.vects["pt_min"] = Point(pt_min)
        self.vects["pt_max"] = Point(pt_max)
        self.inside_points = None
        self.check()

    @property
    def volume(self):
        """Get shape volume"""
        # TODO: The real volume should be computed...
        pt_min = self.vects["pt_min"]()
        pt_max = self.vects["pt_max"]()
        return np.prod(pt_max - p_min)


class Parallelogram(Shape2D):
    def __init__(self, xref, vec1, vec2):
        Shape2D.__init__(self)
        self.vects["xref"] = Point(xref)
        self.vects["vec1"] = Vector(vec1)
        self.vects["vec2"] = Vector(vec2)
        self.check()
        self.bounding_box()

    @property
    def mat(self): return np.vstack((self.vects["vec1"](),
                                     self.vects["vec2"](),))

    def bounding_box(self):
        """Set bounding box for shape"""
        xref = self.vects["xref"]()
        top_pt = xref + self.vects["vec1"]() + self.vects["vec2"]()
        self.vects["pt_min"] = Point([min(xref[0], top_pt[0]), min(xref[1], top_pt[1])])
        self.vects["pt_max"] = Point([max(xref[0], top_pt[0]), max(xref[1], top_pt[1])])

    @property
    def area(self):
        """Get parallelogram area"""
        return abs(np.linalg.det(self.mat))
    volume = area

    def project(self, x):
        """Project vector x on basis formed by parallelogram
        
        x = a * vec1/||vec1|| + b * vec2/||vec2||
        """
        x = np.atleast_2d(x - self.vects["xref"]())
        den = np.linalg.det(self.mat)
        num_a = np.apply_along_axis(
                lambda m: np.linalg.det(np.vstack([m, self.mat[:,1]])), 1, x)
        num_b = np.apply_along_axis(
                lambda m: np.linalg.det(np.vstack([self.mat[:,0], m])), 1, x)
        return num_a / den, num_b / den

    def is_inside(self, x):
        """Determine if x is inside the parallelogram"""
        a, b = self.project(x)
        return (a <= 1) & (a >= 0) & (b <= 1) & (b >= 0)


class Parallelepiped(Shape3D):
    def __init__(self, xref, vec1, vec2, vec3):
        Shape3D.__init__(self)
        self.vects["xref"] = Point(xref)
        self.vects["vec1"] = Vector(vec1)
        self.vects["vec2"] = Vector(vec2)
        self.vects["vec3"] = Vector(vec3)
        self.check()
        self.bounding_box()

    @property
    def mat(self): return np.vstack((self.vects["vec1"](),
                                     self.vects["vec2"](),
                                     self.vects["vec3"](),))

    def bounding_box(self):
        """Set bounding box for shape"""
        xref = self.vects["xref"]()
        top_pt = xref + self.vects["vec1"]() + self.vects["vec2"]() + self.vects["vec3"]()
        self.vects["pt_min"] = Point([min(xref[0], top_pt[0]),
                                      min(xref[1], top_pt[1]),
                                      min(xref[2], top_pt[2])])
        self.vects["pt_max"] = Point([max(xref[0], top_pt[0]),
                                      max(xref[1], top_pt[1]),
                                      max(xref[2], top_pt[2])])

    @property
    def volume(self):
        """Get parallelogram volume"""
        return abs(np.linalg.det(self.mat))

    def project(self, x):
        """Project vector x on basis formed by parallelepiped
        
        x = a * vec1/||vec1|| + b * vec2/||vec2|| + c * vec3/||vec3||
        """
        x = np.atleast_2d(x - self.vects["xref"]())
        den = np.linalg.det(self.mat)
        num_a = np.apply_along_axis(
                lambda m: np.linalg.det(np.vstack([m, self.mat[:,1], self.mat[:,2]])), 1, x)
        num_b = np.apply_along_axis(
                lambda m: np.linalg.det(np.vstack([self.mat[:,0], m, self.mat[:,2]])), 1, x)
        num_c = np.apply_along_axis(
                lambda m: np.linalg.det(np.vstack([self.mat[:,0], self.mat[:,1], m])), 1, x)
        return num_a / den, num_b / den, num_c / den

    def is_inside(self, x):
        """Determine if x is inside the parallelepiped"""
        a, b, c = self.project(x)
        return (a <= 1 ) & (a >= 0) & (b <= 1) & (b >= 0) & (c <= 1) & (c >= 0)


class Disc(Shape2D):
    def __init__(self, center, radius):
        Shape2D.__init__(self)
        self.vects["center"] = Point(center)
        self.radius = radius
        self.check()
        self.bounding_box()

    def bounding_box(self):
        """Set bounding box for shape"""
        self.vects["pt_min"] = self.vects["center"]() - self.radius
        self.vects["pt_max"] = self.vects["center"]() + self.radius

    @property
    def area(self):
        """Get disc area"""
        return np.pi * self.radius**2
    volume = area

    def distance(self, x):
        """Find distance to disc center"""
        return np.linalg.norm(x - self.vects["center"](), axis=-1)

    def is_inside(self, x):
        """Determine if x is inside the disc"""
        return (self.distance(x) <= self.radius)


class Sphere(Shape3D):
    def __init__(self, center, radius):
        Shape3D.__init__(self)
        self.vects["center"] = Point(center)
        self.radius = radius
        self.check()
        self.bounding_box()

    def bounding_box(self):
        """Set bounding box for shape"""
        self.vects["pt_min"] = self.vects["center"]() - self.radius
        self.vects["pt_max"] = self.vects["center"]() + self.radius

    @property
    def volume(self):
        """Get sphere volume"""
        return 4./3. * np.pi * self.radius**3

    def distance(self, x):
        """Find distance to sphere center"""
        return np.linalg.norm(x - self.vects["center"](), axis=-1)

    def is_inside(self, x):
        """Determine if x is inside the sphere"""
        return (self.distance(x) <= self.radius)


class Cylinder(Shape3D):
    def __init__(self, center, radius, vec):
        Shape3D.__init__(self)
        self.vects["xref"] = Point(center)
        self.vects["vec"] = Vector(vec)
        self.radius = radius
        self.check()
        self.bounding_box()

    def bounding_box(self):
        """Set bounding box for shape"""
        def min_max_cyl(axis):
            xref_axis = self.vects["xref"]()[axis]
            vec_axis = self.vects["vec"]()[axis]
            length = np.linalg.norm(self.vects["vec"]())
            proj = self.radius * np.sqrt(1 - vec_axis**2 / length**2)
            four_pts = [xref_axis + proj,
                        xref_axis - proj,
                        xref_axis + vec_axis + proj,
                        xref_axis + vec_axis - proj
                       ]
            return min(four_pts)[0], max(four_pts)[0]
        xmin, xmax = min_max_cyl(0)
        ymin, ymax = min_max_cyl(1)
        zmin, zmax = min_max_cyl(2)
        self.vects["pt_min"] = Point([xmin, ymin, zmin])
        self.vects["pt_max"] = Point([xmax, ymax, zmax])

    @property
    def volume(self):
        """Get cylinder volume"""
        return (np.linalg.norm(self.vects["vec"]()) * np.pi * self.radius**2)[0]

    def project(self, x):
        """Project vector x on basis formed by cylinder
        
        Get cylinder axis coord (from xref) and distance to axis (normalized by ||vec|| and radius)
        x = axis * vec/||vec|| + r * dist_ax/||dist_ax||
        See http://www.flipcode.com/archives/Fast_Point-In-Cylinder_Test.shtml
        """
        x -= self.vects["xref"]()
        dot = np.dot(x, self.vects["vec"]())
        axis_coord = dot / self.vects["vec"].length**2
        r_coord = (np.linalg.norm(x, axis=-1)**2 - axis_coord*dot) / self.radius**2
        return axis_coord, r_coord

    def is_inside(self, x):
        """Determine if x is inside the cylinder"""
        axis, radius = self.project(x)
        return (axis >= 0.) & (axis <= 1.) & (radius <= 1.)
    
    
class rotation():
    """ Quaternions are used to define a rotation on space around any axis
    """
    def __init__(self,origin,direction,theta):
        """
        origin,direction : 3D tuples defining a point and a direction
        theta : angle of rotation in radians
        """
        self.origin = np.asarray(origin)
        self.axis = np.asarray(direction)
        self.theta = theta
        # unit vect as axis
        self.axis /= np.linalg.norm(self.axis)
       
    
    def rotate(self,vect):
        """ vect is a 3D tuple.
        return the rotated vector
        """
        vect0 = np.asarray(vect)
        outvect = np.dot(self._rotation_matrix(),vect0-self.origin)+self.origin
        return outvect
        
        
    def _rotation_matrix(self):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        a = math.cos(self.theta / 2.0)
        b, c, d = -self.axis * math.sin(self.theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

