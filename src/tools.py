#!/usr/bin/env python
"""
Database of flames and subsequent actions

Created Jan 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""

import subprocess

from tempfile import TemporaryFile

def visu(mesh, sol):
    """Call visu script of AVBP (now inside FlameTransfer/bin)"""
    with TemporaryFile() as output:
        process = subprocess.Popen("visu -m {0} -s {1}"
                                   .format(mesh, sol).split(),
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        output, err = process.communicate()
        #while process.wait(): time.sleep(0.01)
        return output
