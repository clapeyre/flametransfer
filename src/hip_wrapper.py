#!/usr/bin/env python
"""
HipWrapper class to handle all executions of hip

Created Mar 2017 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""
import os
import logging
import subprocess

from glob import glob
from time import time, sleep
from tempfile import TemporaryFile

from constants import DEBUG, HIP_START_TIME

class HipWrapper(object):
    """Handle hip executions"""
    def __init__(self, hip_exec):
        self.log = logging.getLogger(__name__)
        self.hip_exec = hip_exec
        self.last_hip_output = ""
        self.last_exec_time = 0.0
        self.execute("", version_check=True)

    def execute(self, script, version_check=False):
        """Execute the current hip script."""
        if version_check:
            arg = "-v"
        else:
            script += "\nqu\n"
            def next_script():
                scripts = glob("script_*.hip")
                if DEBUG:
                    return "script_{:03}.hip".format(len(scripts)+1)
                else:
                    [os.remove(s) for s in scripts]
                    [os.remove(s) for s in glob("script_*.hip.log")]
                    return "script.hip"
            arg = next_script()
            self.log.debug("Executing hip script:")
            [self.log.debug(" > " + line) for line in script.split('\n')]
            with open(arg, 'w') as scr:
                scr.write(script)
        self.last_hip_output = ""
        timestamp = time()
        with TemporaryFile() as output:
            process = subprocess.Popen([self.hip_exec, arg],
                                       stdin=subprocess.PIPE,
                                       stdout=output,
                                       stderr=subprocess.STDOUT)
            while process.poll() is None:
                where = output.tell()
                lines = output.read()
                if not lines:
                    sleep(0.1)
                    output.seek(where)
                else:
                    self.last_hip_output += lines
                if version_check:
                    if time() - timestamp > HIP_START_TIME:
                        process.kill()
                        self.log.error("Hip version was detected as older than "
                                       "17.03, hence unsupported.")
                        self.log.error("Note that if the HIP_START_TIME variable "
                                       "is too low, it can cause false negatives. ")
                        self.log.error("Increase it with `HIP_START_TIME=2.0 "
                                       "flametransfer`")
                        raise AssertionError(
                                "Hip version must be 17.03 at least for "
                                "FlameTransfer. Please upgrade")
            self.last_exec_time = time() - timestamp
            self.log.debug(">>> Hip execution time (s): {}"
                           .format(self.last_exec_time))
            self.last_hip_output += output.read()
            self.log.debug(">>> Hip execution log:")
            for line in self.last_hip_output.split('\n'):
                self.log.debug('\t' + line)
            if process.wait() != 0:
                raise AssertionError("error in hip script. See log")
            if DEBUG:
                with open(path+".log", "w") as out:
                    out.write(self.last_hip_output)
        if not version_check:
            print "\n --- Done executing hip"
