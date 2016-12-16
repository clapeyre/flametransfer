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


class Shape(object):
    @abstractmethod
    def is_inside(self, x):
        pass

    @abstractmethod
    def bounding_box(self):
        pass


class Shape2D(Shape):
    def __init__(self):
        Shape.__init__(self)
        self.ndim = 2
        assert self.area > 1.e-12

    @abstractmethod
    def area(self):
        pass


class Shape3D(Shape):
    def __init__(self):
        Shape.__init__(self)
        self.ndim = 3
        assert self.volume > 1.e-18

    @abstractmethod
    def volume(self):
        pass


class ScatterShape2D(Shape2D):
    def __init__(self, xmin, xmax, ymin, ymax):
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax
        self.inside_points = None
        Shape2D.__init__(self)

    @property
    def area(self):
        """Get shape area"""
        # TODO: The real area should be computed...
        return (self.xmax-self.xmin) * (self.ymax-self.ymin)
    volume = area

    def bounding_box(self):
        """Get bounding box for shape"""
        return {"xmin": self.xmin,
                "ymin": self.ymin,
                "xmax": self.xmax,
                "ymax": self.ymax,
                }


class ScatterShape3D(Shape3D):
    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax):
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.inside_points = None
        Shape3D.__init__(self)

    @property
    def volume(self):
        """Get shape volume"""
        # TODO: The real volume should be computed...
        return (self.xmax-self.xmin) * (self.ymax-self.ymin) * (self.zmax-self.zmin)

    def bounding_box(self):
        """Get bounding box for shape"""
        return {"xmin": self.xmin,
                "ymin": self.ymin,
                "zmin": self.zmin,
                "xmax": self.xmax,
                "ymax": self.ymax,
                "zmax": self.zmax,
                }


class Parallelogram(Shape2D):
    def __init__(self, xref, vec1, vec2):
        self.xref = xref
        self.mat = np.vstack((vec1, vec2))
        Shape2D.__init__(self)

    @property
    def area(self):
        """Get parallelogram area"""
        return abs(np.linalg.det(self.mat))
    volume = area

    def project(self, x):
        """Project vector x on basis formed by parallelogram
        
        x = a * vec1/||vec1|| + b * vec2/||vec2||
        """
        x = np.atleast_2d(x - self.xref)
        den = np.linalg.det(self.mat)
        num_a = np.apply_along_axis(lambda m: np.linalg.det(np.vstack([m, self.mat[:,1]])), 1, x)
        num_b = np.apply_along_axis(lambda m: np.linalg.det(np.vstack([self.mat[:,0], m])), 1, x)
        return num_a / den, num_b / den

    def is_inside(self, x):
        """Determine if x is inside the parallelogram"""
        a, b = self.project(x)
        return (a <= 1) & (a >= 0) & (b <= 1) & (b >= 0)

    def bounding_box(self):
        """Get bounding box for shape"""
        top_pt = self.xref + self.mat[:,0] + self.mat[:,1]
        return {"xmin": min(self.xref[0], top_pt[0]),
                "ymin": min(self.xref[1], top_pt[1]),
                "xmax": max(self.xref[0], top_pt[0]),
                "ymax": max(self.xref[1], top_pt[1])}


class Parallelepiped(Shape3D):
    def __init__(self, xref, vec1, vec2, vec3):
        self.xref = xref
        self.mat = np.vstack((vec1, vec2, vec3))
        Shape3D.__init__(self)

    @property
    def volume(self):
        """Get parallelogram volume"""
        return abs(np.linalg.det(self.mat))

    def project(self, x):
        """Project vector x on basis formed by parallelepiped
        
        x = a * vec1/||vec1|| + b * vec2/||vec2|| + c * vec3/||vec3||
        """
        x = np.atleast_2d(x - self.xref)
        den = np.linalg.det(self.mat)
        num_a = np.apply_along_axis(lambda m: np.linalg.det(np.vstack([m, self.mat[:,1], self.mat[:,2]])), 1, x)
        num_b = np.apply_along_axis(lambda m: np.linalg.det(np.vstack([self.mat[:,0], m, self.mat[:,2]])), 1, x)
        num_c = np.apply_along_axis(lambda m: np.linalg.det(np.vstack([self.mat[:,0], self.mat[:,1], m])), 1, x)
        #num_a = np.linalg.det(np.array((x, self.mat[:,1], self.mat[:,2])).T)
        #num_b = np.linalg.det(np.array((self.mat[:,0], x, self.mat[:,2])).T)
        #num_c = np.linalg.det(np.array((self.mat[:,0], self.mat[:,1]), x).T)
        return num_a / den, num_b / den, num_c / den

    def is_inside(self, x):
        """Determine if x is inside the parallelepiped"""
        a, b, c = self.project(x)
        return (a <= 1 ) & (a >= 0) & (b <= 1) & (b >= 0) & (c <= 1) & (c >= 0)

    def bounding_box(self):
        """Get bounding box for shape"""
        print self.xref
        print self.mat[:,0]
        print self.mat[:,1]
        print self.mat[:,2]
        top_pt = self.xref + self.mat[:,0] + self.mat[:,1] + self.mat[:,2]
        return {"xmin": min(self.xref[0], top_pt[0]),
                "ymin": min(self.xref[1], top_pt[1]),
                "zmin": min(self.xref[2], top_pt[2]),
                "xmax": max(self.xref[0], top_pt[0]),
                "ymax": max(self.xref[1], top_pt[1]),
                "zmax": max(self.xref[2], top_pt[2])}


class Circle(Shape2D):
    def __init__(self, center, radius):
        self.xref = center
        self.radius = radius
        Shape2D.__init__(self)

    @property
    def area(self):
        """Get circle area"""
        return np.pi * self.radius**2
    volume = area

    def distance(self, x):
        """Find distance to circle center"""
        return np.linalg.norm(x - self.xref, axis=-1)

    def is_inside(self, x):
        """Determine if x is inside the circle"""
        return (self.distance(x) <= self.radius)

    def bounding_box(self):
        """Get bounding box for shape"""
        min_point = self.xref - self.radius
        max_point = self.xref + self.radius
        return {"xmin": min_point[0], "xmax": max_point[0],
                "ymin": min_point[1], "ymax": max_point[1]}


class Sphere(Shape3D):
    def __init__(self, center, radius):
        self.xref = center
        self.radius = radius
        Shape3D.__init__(self)

    @property
    def volume(self):
        """Get sphere volume"""
        return 4./3. * np.pi * self.radius**3

    def distance(self, x):
        """Find distance to sphere center"""
        return np.linalg.norm(x - self.xref, axis=-1)

    def is_inside(self, x):
        """Determine if x is inside the circle"""
        return (self.distance(x) <= self.radius)

    def bounding_box(self):
        """Get bounding box for shape"""
        min_point = self.xref - self.radius
        max_point = self.xref + self.radius
        return {"xmin": min_point[0], "xmax": max_point[0],
                "ymin": min_point[1], "ymax": max_point[1],
                "zmin": min_point[2], "zmax": max_point[2]}


class Cylinder(Shape3D):
    def __init__(self, center, radius, vec):
        self.xref = center
        self.radius = radius
        self.vec = vec
        Shape3D.__init__(self)

    @property
    def volume(self):
        """Get cylinder volume"""
        return np.linalg.norm(self.vec) * np.pi * self.radius**2

    def project(self, x):
        """Project vector x on basis formed by cylinder
        
        Get cylinder axis coord and distance to axis (normalized by ||vec|| and radius)
        x = axis * vec/||vec|| + r * dist_ax/||dist_ax||
        See http://www.flipcode.com/archives/Fast_Point-In-Cylinder_Test.shtml
        """
        x -= self.xref
        lengthsqs = np.linalg.norm(self.vec)**2
        dot = np.dot(x, self.vec)
        axis_coord = dot/lengthsqs
        r_coord = (np.linalg.norm(x, axis=-1)**2 - axis_coord*dot) / self.radius**2
        return axis_coord, r_coord

    def is_inside(self, x):
        """Determine if x is inside the cylinder"""
        axis, radius = self.project(x)
        return (axis >= 0.) & (axis <= 1.) & (radius <= 1.)

    def bounding_box(self):
        """Get bounding box for shape"""
        bb1 = Sphere(self.xref, self.radius).bounding_box()
        bb2 = Sphere(self.xref + self.vec, self.radius).bounding_box()
        return {"xmin": min(bb1["xmin"], bb2["xmin"]),
                "ymin": min(bb1["ymin"], bb2["ymin"]),
                "zmin": min(bb1["zmin"], bb2["zmin"]),
                "xmax": max(bb1["xmax"], bb2["xmax"]),
                "ymax": max(bb1["ymax"], bb2["ymax"]),
                "zmax": max(bb1["zmax"], bb2["zmax"])}
    
    
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

