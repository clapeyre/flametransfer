#!/usr/bin/env python
"""
Constants for flametransfer

Created January 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""
__all__ = ['GRID_SIZE',
           'VERSION',
          ]

import os
import sys
import logging
log = logging.getLogger(__name__)

def environ_default(key, default):
    """Get value from environment if defined, otherwise default"""
    env = os.environ.get(key)
    if env is None:
        return default
    else:
        return env

VERSION = "1.0.0-rc.5"
GRID_SIZE = int(environ_default("GRID_SIZE", 32))
DEBUG = environ_default("DEBUG", False)
HIP_START_TIME = float(environ_default("HIP_START_TIME", 1.0))
log.debug("Python executable is: " + sys.executable)

if DEBUG:
    print " !! DEBUG MODE ACTIVATED !!"

try:
    from packaging import version
    CHECK_VERSION = True
except ImportError:
    print " !!WARNING!!: package `packaging` is missing"
    print "              No version checking will be performed"
    CHECK_VERSION = False

class VersionError(Exception):
    """Custom error for version check failure"""
    pass

def compatibility(ver):
    if ver == "1.1B1":
        return "0.3.0-beta"
    return ver

def version_checker(ver):
    log.debug("Checking {} against {}".format(ver, VERSION))
    if not CHECK_VERSION:
        return True
    ver = compatibility(ver)
    ver = version.parse(ver)
    cver = version.parse(VERSION)
    if ver == cver:
        return True
    major, minor, patch = ver.base_version.split('.')
    cmajor, cminor, cpatch = cver.base_version.split('.')
    if major > cmajor:
        raise VersionError(
                "Trying to read object written with future major version. "
                "Please upgrade your software.")
    if major < cmajor:
            print " WARNING : object was written with old discontinued version."
            print "           Compatibility is not guaranteed."
            return False
    elif minor > cminor:
        raise VersionError(
                "Trying to read object written with future version. "
                "Please upgrade your software.")
    elif patch < cpatch:
            print " WARNING : object was written with later (patched) release."
            print "           Your are missing some bugfixes. Update soon."
            return False
    return True
