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

    def help_shell(self):
        print dedent("""\
                Execute a regular shell command
                shell, sh, and ! are equivalent calls""")
    help_sh = help_shell


class FlameTransferCmd(ExitCmd, ShellCmd, cmd.Cmd, object):
    """Command line interpreter for flame transfer"""
    log = logging.getLogger(__name__)
    try:
        hip_exec = os.environ['HIP_EXEC']
    except KeyError:
        print "ERROR: Please define the HIP_EXEC environment variable"
        sys.exit()

    intro = dedent("""\
            Welcome to the FlameTransfer command line
            Type 'help' or '?' for a list of commands, ?about for more on this app""")

    prompt = "ft > "
    flames = []
    current_flame = None

    def help_about(self):
        print dedent("""
                --------------------------------
                Welcome to the FlameTransfer App
                --------------------------------
                
                This app is meant for the import / export and manipulation of
                hdf5 flame files, representing an active flame.  Parameters
                needed for acoustic computations are stored, such as flame
                geometry, reference point and vector, as well as N-tau values.
                Additional data about the flame can also be stored to help
                identify and locate it in the engine.

                Available commands can be listed using 'help' or '?'.  Most
                commands can be abbreviated to the first 2 characters (but not
                'help'). Also, tab completion works, including after a command:

                ft > ge<tab>
                circle         cylinder       parallelogram  sphere 

                Awesome, right?
                """)
    help_ab = help_about

    def emptyline(self):
        """Empty line behavior: do nothing"""
        pass

    generators = "circle cylinder parallelogram sphere".split()
    def do_generate(self, s):
        args = s.split()
        if len(args) != 2:
            print "*** invalid number of arguments"
            return
        if args[1].lower()[:2] == "ci":
            center = input_float(raw_input('Circle center (x y) : '))
            radius = input_float(raw_input('Circle radius (r)   : '))
            self.current_flame = ActiveFlame(args[0], self.hip_exec)
            self.current_flame.define_flame_circle(center, radius)
        if args[1].lower()[:2] == "pa":
            dim = int(raw_input('Number of dimensions (2 or 3) : '))
            if dim == 2:
                xref = input_float(raw_input('Corner point  (x y) : '))
                vec1 = input_float(raw_input('Side vector 1 (x y) : '))
                vec2 = input_float(raw_input('Side vector 2 (x y) : '))
                self.current_flame = ActiveFlame(args[0], self.hip_exec)
                self.current_flame.define_flame_parallelogram(xref, vec1, vec2)
            if dim == 3:
                xref = input_float(raw_input('Corner point  (x y z) : '))
                vec1 = input_float(raw_input('Side vector 1 (x y z) : '))
                vec2 = input_float(raw_input('Side vector 2 (x y z) : '))
                vec3 = input_float(raw_input('Side vector 3 (x y z) : '))
                self.current_flame = ActiveFlame(args[0], self.hip_exec)
                self.current_flame.define_flame_parallelepiped(xref, vec1, vec2, vec3)
            else:
                print "*** nope, it's either 2 or 3"
                return
        if args[1].lower()[:2] == "sp":
            center = input_float(raw_input('Sphere center (x y z) : '))
            radius = input_float(raw_input('Sphere radius (r)     : '))
            self.current_flame = ActiveFlame(args[0], self.hip_exec)
            self.current_flame.define_flame_sphere(center, radius)
        elif args[1].lower()[:2] == "cy":
            center = input_float(raw_input('Cylinder center (x y z) : '))
            radius = input_float(raw_input('Cylinder radius (r)     : '))
            vector = input_float(raw_input('Cylinder vector (x y z) : '))
            self.current_flame = ActiveFlame(args[0], self.hip_exec)
            self.current_flame.define_flame_cylinder(center, radius, vector)
        else:
            print "*** unkown flame shape"
            return
        self.flames.append(self.current_flame)
        self.current_flame.read_meshpoints()
        self.do_shell("rm mesh_{}.*".format(self.current_flame.metas.name))
    do_ge = do_generate

    def help_generate(self):
        print dedent("""\
                Generate an analytical flame shape
                > ge(nerate) <name> <shape>
                |  <name>  : careful, same name flame are overwritten
                |  <shape> : """ + ", ".join(self.generators))
    help_ge = help_generate

    def complete_generate(self, text, line, begidx, endidx):
        if not text:
            completions = self.generators[:]
        else:
            completions = [ f
                            for f in self.generators
                            if f.startswith(text)
                            ]
        return completions
    complete_ge = complete_generate

    writers = "mesh flame ntau".split()
    def do_write(self, s):
        if self.current_flame is None:
            print "*** no flames defined yet"
            return
        if self.current_flame.metas.ref_point is None:
            print "*** please set reference point and vector before writing"
            return
        if s[:2] == "me": # me(sh)
            self.make_mesh()
        elif s[:2] == "fl": # fl(ame)
            self.current_flame.write_h5()
        elif s[:2] == "nt": # nt(au)
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

    def complete_write(self, text, line, begidx, endidx):
        if not text:
            completions = self.writers[:]
        else:
            completions = [ f
                            for f in self.writers
                            if f.startswith(text)
                            ]
        return completions
    complete_wr = complete_write

    listers = "flames refs metas".split()
    def do_list(self, s):
        if self.current_flame is None:
            print "*** no flames defined yet"
            return
        out = open(s.split()[-1], 'w') if (len(s.split()) > 1) else sys.stdout
        if s[:2] == "fl": # fl(ames)
            out.write("     Name       ndim\n")
            for i, f in enumerate(self.flames) :
                star = "*" if f == self.current_flame else " "
                out.write("{0} {1:2} {2:10} {3}\n".format(star, i+1, f.metas.name, f.metas.ndim))
        elif s[:2] == "re": # re(fs)
            out.write("Reference point  : {}\n".format(self.current_flame.metas.ref_point))
            out.write("Reference vector : {}\n".format(self.current_flame.metas.ref_vect))
        elif s[:2] == "me": # me(tas):
            out.write(" Key                  | Value\n")
            out.write(" -------------------------------------------\n")
            for k, v in self.current_flame.metas.__dict__.iteritems():
                out.write(" {0:20} | {1}".format(k, v).replace('\n', ' - ') + '\n')
        else:
            print "*** unknown argument"
            return
    do_li = do_list

    def help_list(self):
        print dedent("""\
                List current memory information
                > li(st) fl(ames)|re(fs) [<path>]
                |   <flames> : view all flames in memory
                |   <refs>   : view reference point and vector of current
                |   <metas>  : view all metas of current
                |   [<path>] : dump output to file <path> instead of stdout""")
    help_li = help_list

    def complete_list(self, text, line, begidx, endidx):
        if not text:
            completions = self.listers[:]
        else:
            completions = [ f
                            for f in self.listers
                            if f.startswith(text)
                            ]
        return completions
    complete_li = complete_list

    setters = "current meta refs ntau".split()
    def do_set(self, s):
        """Change flame meta parameter"""
        if len(s.split()) != 1:
            print "*** invalid number of arguments"
            return
        if self.current_flame is None:
            print "*** no flames defined yet"
            return
        if s[:2] == "cu": # cu(rrent)
            cur = raw_input('Set current flame to : ')
            try:
                nb = int(cur)
                if nb <= len(self.flames):
                    self.current_flame = self.flames[nb-1]
                    return
                else:
                    print "*** there are only {} flames declared".format(len(self.flames))
                    return
            except ValueError:
                for fla in self.flames:
                    if fla.metas.name == cur:
                        self.current_flame = fla
                        return
                print "*** no such flame declared"
                return
        if s[:2] == "me": # me(ta)
            meta = raw_input('Name of meta field to modify : ')
            if not hasattr(self.current_flame.metas, meta):
                print "*** unknown meta"
                return
            typ = raw_input('Type of the field : ')
            val = raw_input('Value             : ')
            if typ == 'int':
                self.current_flame.update_metas({meta: int(val)})
            elif typ == 'float':
                self.current_flame.update_metas({meta: float(val)})
            elif typ == 'array':
                self.current_flame.update_metas({meta: np.array([
                    float(v) for v in val.split()])})
            else:
                print "*** unknown type " + typ
                return
        elif s[:2] == "re": # "re(fs)"
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
            if typ in "n1 crocco cr".split():
                area, p_mean, gamma = input_float(raw_input("Area, Pmean, gamma : "))
                self.current_flame.set_n1_tau(freq, tau, n, area, p_mean, gamma)
            elif typ in "n2 global gl".split():
                self.current_flame.set_n2_tau(freq, tau, n)
            elif typ in "n3 adim ad".split():
                u_bar, q_bar = input_float(raw_input("U_bar, Q_bar : "))
                self.current_flame.set_n3_tau(freq, tau, n, u_bar, q_bar)
        else:
            print "*** unknown settable"
    do_se = do_set

    def help_set(self):
        print dedent("""\
                Set a value or parameter
                > set me(ta) <type> <values(s)>
                |  Set any meta parameter for the flame
                |  <type> = int, float or array
                > set re(fs)
                |  Set the reference point and vector
                > set nt(au)
                |  Set values of N and tau""")
    help_se = help_set

    def complete_set(self, text, line, begidx, endidx):
        if not text:
            completions = self.setters[:]
        else:
            completions = [ f
                            for f in self.setters
                            if f.startswith(text)
                            ]
        return completions
    complete_se = complete_set

    def do_read(self, s):
        args = s.split()
        if len(args) > 2:
            print "*** invalid number of arguments"
            return
        self.current_flame = ActiveFlame('dummy', self.hip_exec)
        if args[-1] == "metas_only":
            self.current_flame.load_metas(path)
        else:
            self.current_flame.load(path)
            self.flames.append(self.current_flame)
    do_re = do_read

    def help_read(self):
        print dedent("""\
                Read hdf5 flame file
                > re(ad) <path> [metas_only]
                | <path>       : absolute or relative path to flame .h5 file
                | [metas_only] : only read the flame's metas. Faster, but incomplete""")
    help_re = help_read


if __name__ == '__main__':
    interpreter = FlameTransferCmd()
    interpreter.cmdloop()
