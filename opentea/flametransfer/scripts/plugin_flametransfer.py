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

        elif action in ["-get_metas-", "-create_flame-", "-duplicate-", "-rename-"]:
            command_exe = self.flametransferexec + " < ./script.ft"
        
        #if command.startswith("-c3sm_auto_loadmesh-"):
        #    scriptFile = os.path.join(self.execDirectory, "./script_mesh.py")
        #    XDR.replace_pattern_in_file(scriptFile, "-c3sm_auto_meshinfo-", os.path.join(self.avbp_toolexec, "meshinfo.e_"+self.avbp_hosttype))
        #    XDR.replace_pattern_in_file(scriptFile, "-c3sm_auto_get_patch_surf-", os.path.join(self.avbp_toolexec, "get_patch_surf.e_"+self.avbp_hosttype))  
        #    command_exe=self.pythonexec   + " ./script_mesh.py"
        #    
        #if command == "-c3sm_auto_preproc_mlpf_exec-":
        #    command_exe = self.avsp_toolexec + "/preproc_mlpf.e_" + self.avsp_hosttype 
        #    
        #if command == "-c3sm_auto_avspinitsol_exec-":
        #        command_exe = self.avsp_toolexec + "/avspinitsol.e_" + self.avsp_hosttype 
        #        
        #if command == "-c3sm_auto_zinn_exec-":
        #    command_exe = self.avsp_toolexec + "/Zinn_nozzle.e_" + self.avsp_hosttype
        #    #print "Warning : the Zinn model exec is not implemented on this plugin"
        #    #return
            
       # if command == "-c3sm_auto_avbp2avsp_exec-":
       #     command_exe = self.avsp_toolexec + "/avbp2avsp_" + self.avsp_hosttype + ".exe"
            
        #if command == "-c3sm_auto_add_activflame-":
        #    global2avsp_cmd = self.avsp_toolexec + "/global2avsp.e_" + self.avsp_hosttype
        #    hipscript_cmd = self.hipexec + " < script_hip"
        #    XDR.replace_pattern_in_file(os.path.join(self.execDirectory, "script_active_flame.py"), "-c3sm_auto_global2avsp_exec-", global2avsp_cmd)
        #    XDR.replace_pattern_in_file(os.path.join(self.execDirectory, "script_active_flame.py"), "-c3sm_auto_hipscript-", hipscript_cmd)
        #    command_exe = self.pythonexec + " ./script_active_flame.py"
        #    
        #if command == "-c3sm_auto_convert_from_avbp-":
        #    avbp2avsp_cmd = self.avsp_toolexec + "/avbp2avsp.e_" + self.avsp_hosttype
        #    hipscript_cmd = self.hipexec + " < script_hip"
        #    XDR.replace_pattern_in_file(os.path.join(self.execDirectory, "script_convert_from_avbp.py"), "-c3sm_auto_avbp2avsp_exec-", avbp2avsp_cmd)
        #    XDR.replace_pattern_in_file(os.path.join(self.execDirectory, "script_convert_from_avbp.py"), "-c3sm_auto_hipscript-", hipscript_cmd)
        #    command_exe = self.pythonexec + " ./script_convert_from_avbp.py"
        #    
        #if command == "-c3sm_auto_duplicate-":
        #    hipscript_cmd = self.hipexec + " < script_hip"
        #    XDR.replace_pattern_in_file(os.path.join(self.execDirectory, "script_duplicate.py"), "-c3sm_auto_hipscript-", hipscript_cmd)
        #    command_exe = self.pythonexec + " ./script_duplicate.py"
        
        else:
            self.log.error("Unknown action " + action)
            return None

        return command_exe 
