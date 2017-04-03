"""
plugin_flametransfer.py

Explicit conversion from generic -c3sm_auto_*- commands to what needs to be done

"""

import logging

class PluginFlametransfer(object):
    """Launch commands"""
    def __init__(self, flametransferexec, hipexec, pythonexec, exec_directory):
        self.flametransferexec = flametransferexec
        self.hipexec = hipexec
        self.pythonexec = pythonexec
        self.exec_directory = exec_directory
        self.log = logging.getLogger(__name__)
        
    def switch_flametransfer_tools(self, action):
        """Switch according to action"""
        if action == "-hip-":
            command_exe = self.hipexec + " ./script.hip"

        elif action in ["-get_metas-", "-create_flame-", "-replicate-", "-rename-"]:
            command_exe = self.flametransferexec + " ./script.ft"
        
        else:
            self.log.error("Unknown action " + action)
            return None

        return command_exe 
