#!/usr/bin/env python
"""
Constants for flametransfer

Created January 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""
__all__ = ['GRID_SIZE',
           'VERSION',
          ]

import os

GRID_SIZE = 32
VERSION = "1.1B1"
DEBUG = True if "DEBUG" in os.environ.keys() else False

#from exceptions import VersionError
class VersionError(Exception):
    pass

def version_checker(version):
    current_ver = version_init(VERSION)
    input_ver = version_init(version)
    if current_ver < input_ver:
        raise VersionError(
                "Trying to read object written with future version. "
                "Please upgrade your software.")
    if current_ver > input_ver:
        raise VersionError(
                "Trying to read object written with outdated version. ")
    if current_ver == input_ver:
        if current_ver.subversion < input_ver.subversion:
            print " WARNING : object was written with later release."
            print "           compatibility should be guaranteed, but please upgrade."
            return False
    return True

def version_init(string):
    if 'B' in string:
        return Beta(string)
    else:
        return Release(string)

class Version(object):
    def __init__(self, string):
        self.string = string

    @property
    def version(self):
        if 'B' in self.string:
            return float(self.string.split('B')[0])
        if 'R' in self.string:
            return float(self.string.split('R')[0])
        return float(self.string)

    def __gt__(self, obj2):
        return self.version > obj2.version

    def __lt__(self, obj2):
        return self.version < obj2.version

    def __ge__(self, obj2):
        return self.version >= obj2.version

    def __le__(self, obj2):
        return self.version <= obj2.version

    def __eq__(self, obj2):
        return self.version == obj2.version

class Release(Version):
    @property
    def subversion(self):
        if 'R' in self.string:
            return int(self.string.split('R')[1])
        return 0

class Beta(Version):
    @property
    def version(self):
        version, subversion = self.string.split('B')
        return float(version) + 0.001*float(subversion)

    @property
    def subversion(self): return 0
