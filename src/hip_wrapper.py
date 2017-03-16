#!/usr/bin/env python
"""
HipWrapper class to handle all executions of hip

Created Mar 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""
import os
import time
import logging
import subprocess

from glob import glob
from tempfile import TemporaryFile

from constants import DEBUG

class HipWrapper(object):
    """Handle hip executions"""
    def __init__(self, hip_exec):
        self.log = logging.getLogger(__name__)
        self.hip_exec = hip_exec
        self.last_hip_output = ""

    def execute(self, script):
        """Execute the current hip script."""
        script += "\nqu\n"
        def next_script():
            scripts = glob("script_*.hip")
            if DEBUG:
                return "script_{:03}.hip".format(len(scripts)+1)
            else:
                [os.remove(s) for s in scripts]
                [os.remove(s+".log") for s in scripts]
                return "script.hip"
        path = next_script()
        self.last_hip_output = ""
        self.log.debug("Executing hip script:")
        [self.log.debug(" > " + line) for line in script.split('\n')]
        with open(path, 'w') as scr:
            scr.write(script)
        with TemporaryFile() as output:
            process = subprocess.Popen([self.hip_exec, path],
                                       stdin=subprocess.PIPE,
                                       stdout=output,
                                       stderr=subprocess.STDOUT)
            while process.poll() is None:
                where = output.tell()
                lines = output.read()
                if not lines:
                    time.sleep(0.1)
                    output.seek(where)
                else:
                    self.last_hip_output += lines
            self.last_hip_output += output.read()
            self.log.debug("Hip execution log")
            for line in self.last_hip_output.split('\n'):
                self.log.debug(line)
            if process.wait() != 0:
                raise AssertionError("error in hip script. See log")
            if DEBUG:
                with open(path+".log", "w") as out:
                    out.write(self.last_hip_output)
        print "\n --- Done executing hip"
