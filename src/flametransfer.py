#!/usr/bin/env python
"""
FlameTransfer main command line utility

Created December 2016 by Corentin Lapeyre (lapeyre@cerfacs.fr)
"""

import os
import logging
import sys
import cmd
import copy
import numpy as np

from os.path import isfile
from textwrap import dedent
from StringIO import StringIO

from activeflame import ActiveFlame
from geometry import NormalVector, Vector, Point
from constants import VERSION

logger = logging.getLogger()
logger.setLevel("DEBUG")
file_format = logging.Formatter("%(asctime)s %(name)s %(levelname)s  %(message)s")
file_handler = logging.FileHandler('flametransfer.log')
file_handler.setFormatter(file_format)
logger.addHandler(file_handler)

def input_floats(inp):
    """Convert raw input to float array"""
    return [float(f) for f in inp.split()]

def input_array(inp):
    """Convert raw input to numpy array"""
    return np.array(input_floats(inp))

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
        print "\n --- Exiting FlameTransfer. Bye!"
        sys.exit(0)
        return True
    do_EOF = do_exit
    do_quit = do_exit


class ShellCmd(cmd.Cmd, object):
    def do_shell(self, s):
        """Execute a regular shell command"""
        os.system(s)
    do_sh = do_shell

    def help_shell(self):
        print dedent("""\
                Execute a regular shell command. 
                Useful for e.g. 'shell ls' (to see what has been written).
                Note : '!ls' is equivalent to 'shell ls'.
                Warning : Your .bashrc file is *not* sourced.""")
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
            Welcome to the FlameTransfer V.{} command line
            Type help or ? for a list of commands,
                 ?about for more on this app""").format(VERSION)

    prompt = "ft > "
    flames = []
    cur_fl = None

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

                Available commands can be listed using 'help' or '?'.  All
                commands can be abbreviated to the shortest character list that
                does not lead to ambiguity (2 is often good, 3 should be very
                safe). Everything is parsed based on whitespaces, so a vector
                is simply:
                
                Please input vector : 1 0 0
                
                This implies that 'my new flame' will not work as a flame name.
                Also, tab completion works, including after a command:

                ft > gen<tab>
                ft > generate <tab>
                avbp_field     circle         cylinder       parallelogram  sphere 

                Awesome, right?
                """)

    def default(self, line):
        """Override this command from cmd.Cmd to accept shorcuts"""
        cmd, arg, _ = self.parseline(line)
        func = [getattr(self, n) for n in self.get_names() if n.startswith('do_' + cmd)]
        if len(func) == 0:
            self.stdout.write('*** Unknown syntax: %s\n'%line)
            return
        elif len(func) > 1:
            print '*** {} is a shorcut to several commands'.format(cmd)
            print '    Please give more charaters for disambiguation'
            return
        else:
            func[0](arg)

    def do_help(self, arg):
        """Wrapper for cmd.Cmd.do_help to accept shortcuts"""
        if arg:
            helper = [n[5:] for n in self.get_names() if n.startswith('help_' + arg)]
            if len(helper) == 0:
                self.stdout.write('*** Unknown command: %s\n'%arg)
                return
            elif len(helper) > 1:
                self.stdout.write('*** {} is a shorcut to several commands'.format(cmd))
                self.stdout.write('    Please give more charaters for disambiguation')
                return
            else:
                arg = helper[0]
        cmd.Cmd.do_help(self, arg) 

    def emptyline(self):
        """Empty line behavior: do nothing"""
        pass

    generators = "avbp_field disc cylinder parallelogram sphere".split()
    def do_generate(self, s):
        args = s.split()
        if len(args) != 2:
            print "*** invalid number of arguments"
            return
        typ, name = args
        typ = typ[0].lower()
        if typ == "d": # disc
            center = input_floats(raw_input('Circle center (x y) : '))
            radius = input_floats(raw_input('Circle radius (r)   : '))[0]
            self.cur_fl = ActiveFlame(name, self.hip_exec)
            self.cur_fl.define_flame_disc(center, radius)
        elif typ == "p": # parallelogram
            dim = int(raw_input('Number of dimensions (2 or 3) : '))
            if dim == 2:
                xref = input_floats(raw_input('Corner point  (x y) : '))
                vec1 = input_floats(raw_input('Side vector 1 (x y) : '))
                vec2 = input_floats(raw_input('Side vector 2 (x y) : '))
                self.cur_fl = ActiveFlame(name, self.hip_exec)
                self.cur_fl.define_flame_parallelogram(xref, vec1, vec2)
            elif dim == 3:
                xref = input_floats(raw_input('Corner point  (x y z) : '))
                vec1 = input_floats(raw_input('Side vector 1 (x y z) : '))
                vec2 = input_floats(raw_input('Side vector 2 (x y z) : '))
                vec3 = input_floats(raw_input('Side vector 3 (x y z) : '))
                self.cur_fl = ActiveFlame(name, self.hip_exec)
                self.cur_fl.define_flame_parallelepiped(xref, vec1, vec2, vec3)
            else:
                print "*** nope, it's either 2 or 3"
                return
        elif typ == "s": # sphere
            center = input_floats(raw_input('Sphere center (x y z) : '))
            radius = input_floats(raw_input('Sphere radius (r)     : '))
            self.cur_fl = ActiveFlame(name, self.hip_exec)
            self.cur_fl.define_flame_sphere(center, radius)
        elif typ == "c": # cylinder
            center = input_floats(raw_input('Cylinder center (x y z) : '))
            radius = input_floats(raw_input('Cylinder radius (r)     : '))
            vector = input_floats(raw_input('Cylinder vector (x y z) : '))
            self.cur_fl = ActiveFlame(name, self.hip_exec)
            self.cur_fl.define_flame_cylinder(center, radius, vector)
        elif typ == "a": # avbp field
            avbp_mesh = raw_input('AVBP mesh : ')
            if not isfile(avbp_mesh):
                print "*** not such file"
                return
            avbp_sol = raw_input('AVBP sol : ')
            if not isfile(avbp_sol):
                print "*** not such file"
                return
            avbp_field = raw_input('Scalar field : ')
            avbp_thresh = raw_input('Threshold value : ')
            self.cur_fl = ActiveFlame(name, self.hip_exec)
            self.cur_fl.define_threshold_flame(
                    avbp_mesh, avbp_sol, avbp_field, avbp_thresh)
        else:
            print "*** unkown flame shape"
            return
        self.flames.append(self.cur_fl)
        self.cur_fl.read_meshpoints()
        #self.do_shell("rm mesh_{}.*".format(self.cur_fl.metas.name))

    def help_generate(self):
        print dedent("""\
                Generate an analytical flame shape. All lenths are in [m]
                > ge(nerate) <shape> <name> 
                |  <shape> : one of {}
                |  <name>  : careful, same name flame is overwritten""".format(
                  ", ".join(self.generators)))

    def complete_generate(self, text, line, begidx, endidx):
        return [ f for f in self.generators if f.startswith(text) ]

    writers = "mesh flame full ntau all".split()
    def do_write(self, s):
        if self.cur_fl is None:
            print "*** no flames defined yet"
            return
        if self.cur_fl.metas.pt_ref is None:
            print "*** please set reference point and vector before writing"
            return
        if s[:2] == "me": # me(sh)
            self.cur_fl.make_mesh()
        elif s[:2] == "fl": # fl(ame)
            self.cur_fl.write_h5()
        elif s[:2] == "fu": # fu(ll)
            self.cur_fl.write_h5(with_mesh=True)
        elif s[:2] == "nt": # nt(au)
            if self.cur_fl.metas.n2_tau is None:
                print "*** please set N-tau before writing"
                return
            self.cur_fl.write_n_tau()
        elif s[:2] == "al": # al(l)
            self.cur_fl.write_h5(with_mesh=True)
            if self.cur_fl.metas.n2_tau is None:
                print "*** please set N-tau before writing"
                return
            self.cur_fl.write_n_tau()
        else:
            print "*** unknown option"
            return

    def help_write(self):
        print dedent("""\
                Write output file
                > wr(ite) me(sh)|fl(ame)|nt(au)
                |  mesh  : write hdf5 mesh file
                |  flame : write hdf5 flame file
                |  full  : write hdf5 mesh, flame and xmf file
                |  ntau  : write ascii n-tau file
                |  all   : write all of the above""")

    def complete_write(self, text, line, begidx, endidx):
        return [ f for f in self.writers if f.startswith(text) ]

    listers = "flames refs metas json dump".split()
    def do_list(self, s):
        if self.cur_fl is None:
            print "*** no flames defined yet"
            return
        out = open(s.split()[-1], 'w') if (len(s.split()) > 1) else sys.stdout
        if out == sys.stdout: print
        arg = s[0].lower()
        if arg == "f": # flames
            width = max([len(f.metas.name) for f in self.flames])
            out.write("Cu   Name" + " "*width + "ndim  Volume     PtRef           VecRef\n")
            out.write("-" * (48 + width) + "\n")
            for i, f in enumerate(self.flames) :
                star = "*" if f == self.cur_fl else " "
                out.write("{0} {1:2} {2:{width}}    {3}   {4:.2E}   {5: <15} {6}\n".format(
                    star, i+1, f.metas.name, f.metas.ndim, f.shape.volume,
                    f.metas.pt_ref, f.metas.vec_ref, width=width+2))
        elif arg == "r": # refs
            out.write("Reference point  : {}\n".format(self.cur_fl.metas.pt_ref))
            out.write("Reference vector : {}\n".format(self.cur_fl.metas.vec_ref))
        elif arg == "m": # metas
            out.write(" Key                  | Value\n")
            out.write(" -------------------------------------------\n")
            for k, v in sorted(self.cur_fl.metas.static.items()):
                out.write(" {0:20} | {1}".format(k, v).replace('\n', ' ') + '\n')
            for k, v in sorted(self.cur_fl.metas.vects.items()):
                out.write(" {0:20} | {1}".format(
                    k, v() if v is not None else None).replace('\n', ' ') + '\n')
        elif arg == "j": # json formatted output
            try:
                with open("flame_{}.metas.json".format(
                              self.cur_fl.metas.name),
                          'w') as f:
                    self.cur_fl.metas.write_json(f)
            except ImportError:
                print "*** could not import python module json"
                return
        elif arg == "d": # dump
            if out == sys.stdout:
                print "*** <key> argument mandatory for dump"
                return
            key = s.split()[-1]
            if key not in self.cur_fl.metas.__dict__.keys():
                print "*** unknown key " + key
                return
            out.write("{}".format(
                self.cur_fl.metas.__dict__[key]).translate(None, '[],'))
        else:
            print "*** unknown argument"
            return
        if out != sys.stdout:
            out.close()

    def help_list(self):
        print dedent("""\
                List current memory information
                > li(st) fl(ames)|re(fs)|me(tas) [<path>]
                |   flames   : view all flames in memory
                |   refs     : view reference point and vector of current
                |   metas    : view all metas of current
                |   [<path>] : dump output to file <path> instead of stdout
                > li(st) js(on)
                |   No arguments, a flame_{name}.metas.json file is written with all metas
                > li(st) du(mp) <key>
                |   <key>    : input a key, get the resulting meta value in a '<key>' file""")

    def complete_list(self, text, line, begidx, endidx):
        return [ f for f in self.listers if f.startswith(text) ]

    setters = "current meta refs ntau".split()
    def do_set(self, s):
        """Change flame meta parameter"""
        if len(s.split()) != 1:
            print "*** invalid number of arguments"
            return
        if self.cur_fl is None:
            print "*** no flames defined yet"
            return
        switch = s[0].lower()
        if switch == "c": # current
            cur = raw_input('Set current flame to : ')
            try:
                nb = int(cur)
                if nb <= len(self.flames):
                    self.cur_fl = self.flames[nb-1]
                    return
                else:
                    print "*** there are only {} flames declared".format(len(self.flames))
                    return
            except ValueError:
                for fla in self.flames:
                    if fla.metas.name == cur:
                        self.cur_fl = fla
                        return
                print "*** no such flame declared"
                return
        elif switch == "s": # static meta
            meta = raw_input('Name of meta field to set : ')
            typ = raw_input('Type of the field : ')[0].lower()
            typ_dict = {'i': int,
                        'f': float,
                        'a': input_floats,
                        's': lambda s: s,} # there is no 'identity' func in python!
            if typ not in typ_dict.keys():
                print "*** unknown type " + typ
                return
            val = typ_dict[typ](raw_input('Value             : '))
            setattr(self.cur_fl.metas, meta, val)
        elif switch == "v": #vector meta
            meta = raw_input('Name of meta field to set : ')
            typ = raw_input('Type of the field : ')[0].lower()
            typ_dict = {'n': NormalVector,
                        'v': Vector,
                        'p': Point,}
            if typ not in typ_dict.keys():
                print "*** unknown type " + typ
                return
            self.cur_fl.metas.set_vect(
                    meta,
                    typ_dict[typ](input_floats(raw_input('Value             : '))))
        elif switch == "r": # refs
            metas = {}
            self.cur_fl.metas.pt_ref = Point(input_floats(raw_input('Reference point  (x y z) : ')))
            self.cur_fl.metas.vec_ref = Vector(input_floats(raw_input('Reference vector (x y z) : ')))
        elif switch == "n": # ntau
            typ = raw_input("Type of N : ").lower()
            if typ not in "n1 n2 n3 crocco global adim".split():
                print "*** unknown N type"
                return
            path = raw_input("Path to N-tau file : ")
            def exp_wrapper(path):
                with open(path, 'r') as f:
                    return StringIO(f.read().lower().replace('d', 'e'))
            freq, n, tau = np.loadtxt(exp_wrapper(path)).T
            typ = typ[:2].lower()
            if typ in "n1 cr".split():
                area, p_mean, gamma = input_floats(raw_input("Area, Pmean, gamma : "))
                self.cur_fl.set_n1_tau(freq, tau, n, area, p_mean, gamma)
            elif typ in "n2 gl".split():
                self.cur_fl.set_n2_tau(freq, tau, n)
            elif typ in "n3 ad".split():
                u_bar, q_bar = input_floats(raw_input("U_bar, Q_bar : "))
                self.cur_fl.set_n3_tau(freq, tau, n, u_bar, q_bar)
            else:
                print "*** unknown N type"
                return
        else:
            print "*** unknown settable"

    def help_set(self):
        print dedent("""\
                Set a value or parameter for the current flame
                > set current
                |  Set current flame to <ref>
                |  <ref> = either a flame name or number
                |  See ft> list flames to view which is current
                > set static
                |  Set any meta parameter for the flame
                |  types: string, int, float or array
                > set vector
                |  Set a point or vector meta parameter for the flame
                |  types: normal_vector, vector, point
                > set refs
                |  Set the reference point and vector
                > set ntau
                |  Set values of N and tau
                |  Accepted N types: {n1|crocco, n2|global, n3|adim}
                |  See ?N_tau for more""")

    def complete_set(self, text, line, begidx, endidx):
        return [ f for f in self.setters if f.startswith(text) ]

    def help_N_tau(self):
        print dedent("""
                ----------
                N-tau data
                ----------
                
                In AVSP 5.6 and above, N-tau data is frequency dependant.
                Input files are column tab/space separated ascii files, and
                each line contains:
                 
                freq N tau

                '#' is the comment character, so feel free to add comments to
                identify your data.  The number of tabs or spaces used does not
                matter.

                AVSP then interpolates linearly in this table.  First and
                last values are kept constant:
                |         .___. 
                |       ./     \.___
                |  ___./  
                Example: if this is my data file :
                
                # My great data
                100 1.1 0.001
                500 0.9 0.003  # Not sure about the delay
                
                then AVSP will compute for the following :
                | 50 Hz: N=1.1, tau=0.001
                |300 Hz: N=1.0, tau=0.002
                |900 Hz: N=0.9, tau=0.003
                """)

    def do_read(self, s):
        args = s.split()
        if len(args) > 2 or len(args) < 1:
            print "*** invalid number of arguments"
            return
        self.cur_fl = ActiveFlame('tmp', self.hip_exec)
        if not isfile(args[0]):
            print "*** no such file"
            return
        if args[-1] == "metas_only":
            self.cur_fl.load_metas(args[0])
        else:
            self.cur_fl.load(args[0])
            self.flames.append(self.cur_fl)

    def help_read(self):
        print dedent("""\
                Read hdf5 flame file
                > re(ad) <path> [metas_only]
                |  <path>       : absolute or relative path to flame .h5 file
                |  [metas_only] : only read the flame's metas. Faster, but incomplete""")

    def do_copy(self, s):
        args = s.split()
        if len(args) != 2:
            print "*** invalid number of arguments"
            return
        if self.cur_fl is None:
            print "*** no flames defined yet"
            return
        src, dest = args
        if dest in [f.metas.name for f in self.flames]:
            print "*** {} is already a declared flame".format(dest)
            return
        try:
            nb = int(src)
            if nb <= len(self.flames):
                self.flames.append(copy.deepcopy(self.flames[nb-1]))
                self.cur_fl = self.flames[-1]
                return
            else:
                print "*** there are only {} flames declared".format(len(self.flames))
                return
        except ValueError:
            for fla in self.flames:
                if fla.metas.name == src:
                    self.flames.append(copy.deepcopy(fla))
                    self.cur_fl = self.flames[-1]
                    return
            print "*** no such flame declared"
            return

    def help_copy(self):
        print dedent("""\
                Copy existing flame
                > co(py) <name_src> <name_dest>
                |  <name_src>  : name or number of flame to copy
                |  <name_dest> : name of destination flame""")

    transforms = "translate scale rotate".split()
    def do_transform(self, s):
        """Appy transformation to current flame"""
        if self.cur_fl is None:
            print "*** no flames defined yet"
            return
        switch = s[0].lower()
        if switch == "t":
            vect = input_floats(raw_input('Translation vector : '))
            self.cur_fl.transform(switch, vect)
        elif switch == "s":
            vect = input_floats(raw_input('Scale vector : '))
            self.cur_fl.transform(switch, vect)
        elif switch == "r":
            axis = raw_input("Rotation axis : ")
            angle = raw_input("Rotation angle (degrees) : ")
            self.cur_fl.transform(switch, (axis, float(angle)))
        else:
            print "*** Accepted transforms : " + " ".join(self.transforms)
            return

    def help_transform(self):
        print dedent("""\
                Apply transformation to current flame
                > transform translate|scale|rotate""")

    def complete_transform(self, text, line, begidx, endidx):
        return [ f for f in self.transforms if f.startswith(text) ]

    exports = "AVSP5.6".split()
    def do_export(self, s):
        """Export current flame"""
        if self.cur_fl is None:
            print "*** no flames defined yet"
            return
        if s.lower() == "avsp5.6":
            avsp_mesh = raw_input('AVSP mesh file : ')
            avsp_sol = raw_input('AVSP sol file : ')
            self.cur_fl.export_avsp(avsp_mesh, avsp_sol)
        else:
            print "*** unknown export format"

    def help_export(self):
        print dedent("""\
                Export current flame to AVSP solution
                > exp(ort) AVSP5.6
                |  variable frequency n-tau but single flame AVSP version""")

    def complete_export(self, text, line, begidx, endidx):
        return [ f for f in self.exports if f.startswith(text) ]

if __name__ == '__main__':
    if len(sys.argv) > 1:
        print "*** to execute a script, please feed to stdin:"
        print "    $ flametransfer < " + sys.argv[1]
        sys.exit()
    interpreter = FlameTransferCmd()
    interpreter.cmdloop()
