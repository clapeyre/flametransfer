#!/usr/bin/env python
"""
FlameTransfer main command line utility

Created December 2016 by COOP team
"""

import os
import logging
import sys
import cmd
import numpy as np

from textwrap import dedent

from activeflame import ActiveFlame

logger = logging.getLogger()
logger.setLevel("DEBUG")
file_format = logging.Formatter("%(asctime)s %(name)s %(levelname)s  %(message)s")
file_handler = logging.FileHandler('flametransfer.log')
file_handler.setFormatter(file_format)
logger.addHandler(file_handler)

def input_float(inp):
    """Convert raw input to numpy array"""
    return np.array([float(f) for f in inp.split()])

class ExitCmd(cmd.Cmd, object):
    def can_exit(self):
        """This can be changed if exit must be protected"""
	return True

    def onecmd(self, line):
	r = super (ExitCmd, self).onecmd(line)
	if r and (self.can_exit() or
	   raw_input('exit anyway ? (yes/no):')=='yes'):
	     return True
	return False

    def do_exit(self, s):
        """Exit the interpreter."""
        """You can also use the Ctrl-D shortcut."""
	return True
    do_EOF = do_exit
    do_ex = do_exit
    do_quit = do_exit
    do_qu = do_exit


class ShellCmd(cmd.Cmd, object):
    def do_shell(self, s):
        """Execute a regular shell command"""
        os.system(s)
    do_sh = do_shell


class FlameTransferCmd(ExitCmd, ShellCmd, cmd.Cmd, object):
    """Command line interpreter for flame transfer"""
    log = logging.getLogger(__name__)
    try:
        hip_exec = os.environ['HIP_EXEC']
    except KeyError:
        print "ERROR: Please define the HIP_EXEC environment variable"
        sys.exit()

    intro = dedent("""\
            Welcome to the flametransfer command line
            Type 'help' for a list of commands""")

    prompt = "ft > "
    flames = []
    current_flame = None

    def exec_hip(self):
        """Execute the current hip script."""
        if not self.current_flame.do_hip:
            print "*** flame is not ready to execute hip"
            return
        self.do_shell(self.hip_exec + ' < script.hip; echo')
        self.current_flame.do_hip = False

    def help_introduction(self):
        print 'introduction'
        print 'a good place for a tutorial'

    def emptyline(self):
        """Empty line behavior: do nothing"""
        pass

    def do_add(self, s):
        l = s.split()
        if len(l)!=2:
           print "*** invalid number of arguments"
           return
        try:
           l = [int(i) for i in l]
        except ValueError:
           print "*** arguments should be numbers"
           return
        print l[0]+l[1]

    def do_generate(self, s):
        args = s.split()
        if len(args) != 2:
            print "*** invalid number of arguments"
            return
        if args[1].lower()[:2] == "ci":
            center = input_float(raw_input('Circle center (x y) : '))
            radius = input_float(raw_input('Circle radius (r)   : '))
            self.current_flame = ActiveFlame(args[0])
            self.current_flame.define_flame_circle(center, radius)
        if args[1].lower()[:2] == "pa":
            xref = input_float(raw_input('Corner point  (x y) : '))
            vec1 = input_float(raw_input('Side vector 1 (x y) : '))
            vec2 = input_float(raw_input('Side vector 2 (x y) : '))
            self.current_flame = ActiveFlame(args[0])
            self.current_flame.define_flame_parallelogram(xref, vec1, vec2)
        if args[1].lower()[:2] == "sp":
            center = input_float(raw_input('Sphere center (x y z) : '))
            radius = input_float(raw_input('Sphere radius (r)     : '))
            self.current_flame = ActiveFlame(args[0])
            self.current_flame.define_flame_sphere(center, radius)
        elif args[1].lower()[:2] == "cy":
            center = input_float(raw_input('Cylinder center (x y z) : '))
            radius = input_float(raw_input('Cylinder radius (r)     : '))
            vector = input_float(raw_input('Cylinder vector (x y z) : '))
            self.current_flame = ActiveFlame(args[0])
            self.current_flame.define_flame_cylinder(center, radius, vector)
        else:
            print "*** unkown flame shape"
            return
        self.flames.append(self.current_flame)
        self.exec_hip()
        self.current_flame.read_meshpoints()
        self.do_shell("rm " + self.current_flame.mesh_file)
    do_ge = do_generate

    def do_write(self, s):
        if self.current_flame is None:
            print "*** no flames defined yet"
            return
        if self.current_flame.metas.ref_point is None:
            print "*** please set reference point and vector before writing"
            return
        if s[:2] == "me":
            self.make_mesh_script()
            self.exec_hip()
        elif s[:2] == "fl":
            self.current_flame.write_h5()
        elif s[:2] == "nt":
            if self.current_flame.metas.n2_tau is None:
                print "*** please set N-tau before writing"
                return
            self.current_flame.write_n_tau()
        else:
            print "*** unknown option"
            return
    do_wr = do_write

    def help_write(self):
        print dedent("""\
                Write output file
                > wr(ite) me(sh)|fl(ame)|nt(au)
                |  mesh  : write hdf5 mesh file
                |  flame : write hdf5 flame file
                |  ntau  : write ascii n-tau file""")
    help_wr = help_write

    def help_generate(self):
        print dedent("""\
                Generate an analytical flame shape
                > ge(nerate) <name> <shape>
                |  <name>  : careful, same name flame are overwritten
                |  <shape> : sphere, cylinder""")
    help_ge = help_generate

    def do_list(self, s):
        if self.current_flame is None:
            print "*** no flames defined yet"
            return
        if s[:2] == "fl":
            print "     Name       ndim"
            for i, f in enumerate(self.flames) :
                star = "*" if f == self.current_flame else " "
                print "{0} {1:2} {2:10} {3}".format(star, i+1, f.metas.name, f.metas.ndim)
        elif s[:2] == "re":
            print "Reference point  :", self.current_flame.metas.ref_point
            print "Reference vector :", self.current_flame.metas.ref_vect
        else:
            print "*** unknown argument"
    do_li = do_list

    def help_list(self):
        print dedent("""\
                List current memory information
                > li(st) fl(ames)|re(fs)
                |   <flames> : view all flames in memory
                |   <refs>   : view reference point and vector""")
    help_li = help_list

    def do_set(self, s):
        """Change flame meta parameter"""
        if len(s.split()) != 1:
            print "*** invalid number of arguments"
            return
        if self.current_flame is None:
            print "*** no flames defined yet"
            return
        if s[:2] == "me":
            meta = raw_input('Name of meta field to modify : ')
            if not hasattr(self.current_flame.metas, meta):
                print "*** unknown meta"
                return
            typ = raw_input('Type of the field : ')
            if typ == 'int':
                self.current_flame.update_metas({args[1]: int(args[1])})
            elif typ == 'float':
                self.current_flame.update_metas({args[1]: float(args[1])})
            elif typ == 'array':
                self.current_flame.update_metas({args[1]: np.array([
                    float(v) for v in args[3:]])})
            else:
                print "*** unknown type " + args[2]
                return
        elif s[:2] == "re":
            metas = {}
            metas["ref_point"] = input_float(raw_input('Reference point  (x y z) : '))
            metas["ref_vect"] = input_float(raw_input('Reference vector (x y z) : '))
            self.current_flame.update_metas(**metas)
        elif s[:2] == "nt":
            typ = raw_input("Type of N : ").lower()
            if typ not in "n1 n2 n3 crocco global adim".split():
                print "*** unknown N type"
                return
            path = raw_input("Path to N-tau file : ")
            freq, n, tau = np.loadtxt(path)
            if typ in "n1 crocco".split():
                area, p_mean, gamma = input_float(raw_input("Area, Pmean, gamma : "))
                self.current_flame.set_n1_tau(freq, tau, n, area, p_mean, gamma)
            elif typ in "n2 global".split():
                self.current_flame.set_n2_tau(freq, tau, n)
            elif typ in "n3 adim".split():
                u_bar, q_bar = input_float(raw_input("U_bar, Q_bar : "))
                self.current_flame.set_n3_tau(freq, tau, n, u_bar, q_bar)
        else:
            print "*** unknown settable"
    do_se = do_set

    def help_set(self):
        print dedent("""\
                Set a value or parameter
                > set meta <type> <values(s)>
                |  Set any meta parameter for the flame
                |  <type> = int, float or array
                > set refs
                |  Set the reference point and vector""")
    help_se = help_set

if __name__ == '__main__':
    interpreter = FlameTransferCmd()
    interpreter.cmdloop()
