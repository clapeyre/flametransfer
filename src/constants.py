#!/usr/bin/env python
"""
Constants for flametransfer

Created January 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""
__all__ = ['GRID_SIZE',
           'VERSION',
          ]

import os
import semver

def environ_default(key, default):
    """Get value from environment if defined, otherwise default"""
    env = os.environ.get(key)
    if env is None:
        return default
    else:
        return env

VERSION = "1.1.0-alpha"
GRID_SIZE = int(environ_default("GRID_SIZE", 32))
DEBUG = environ_default("DEBUG", False)
HIP_START_TIME = float(environ_default("HIP_START_TIME", 1.0))

class VersionError(Exception):
    """Custom error for version check failure"""
    pass

def version_checker(version):
    """Check input version vs current using semantic versioning"""
    if version == "1.1B1":
        version = "1.0.1-beta"
    current_ver = semver.parse_version_info(VERSION)
    input_ver = semver.parse_version_info(version)
    if current_ver < input_ver:
        raise VersionError(
            "Trying to read object written with future version. "
            "Please upgrade your flametransfer software.")
    if current_ver.major > input_ver.major:
        raise VersionError(
            "Flametransfer {0}.X.X cannot read old data from {1}"
            .format(current_ver.major, version))
    return True
