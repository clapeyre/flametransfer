#!/usr/bin/env python
"""
Database of flames and subsequent actions

Created Jan 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""

from visu import main as visu_main


def visu(mesh, sol):
    """Call visu script of AVBP (now inside FlameTransfer/bin)"""
    visu_main(mesh, sol, None, None, None)
