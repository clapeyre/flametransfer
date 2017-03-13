#!/usr/bin/env python
"""
FlameTransfer main command line utility

USAGE:
    There are 2 ways to interact with this script:
    - as an interactive command line. Type `flametransfer` and make all the
      mistakes that you want, errors will be caught and explained.
    - for batch processing. Write a file `cmd` containing commands as you
      would enter them on the interactive command line. Errors will be raised
      and the script interrupted.

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
from flamedb import FlameDB

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

def is_number(inp):
    """Check if string is a number"""
    try:
        float(inp)
        return True
    except ValueError:
        return False


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


class SmartCmd(cmd.Cmd, object):
    """Good featured command line
    
    - function / help shortcuts (as short as disambiguation permits)
    - catch ^C
    - catch assertion errors

    Note: behavior depends on self.use_rawinput:
    > if True: interactive session expected. Catch errors
    > if False: batch session, raise exceptions and interrupt
    """
    def cmdloop_with_keyboard_interrupt(self):
        doQuit = False
        while doQuit != True:
            try:
                self.cmdloop()
                doQuit = True
            except KeyboardInterrupt:
                self.intro = None
                sys.stdout.write('\n')

    def default(self, line):
        """Override this command from cmd.Cmd to accept shorcuts"""
        cmd, arg, _ = self.parseline(line)
        func = [getattr(self, n) for n in self.get_names() if n.startswith('do_' + cmd)]
        if len(func) == 0:
            self.stdout.write('*** unknown syntax: %s\n'%line)
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
                self.stdout.write('*** unknown command: %s\n'%arg)
                return
            elif len(helper) > 1:
                self.stdout.write('*** {} is a shorcut to several commands'.format(cmd))
                self.stdout.write('    Please give more charaters for disambiguation')
                return
            else:
                arg = helper[0]
        cmd.Cmd.do_help(self, arg) 

    def onecmd(self, line):
        """Wrapper for cmd.Cmd.onecmd to catch assertion errors"""
        try:
            cmd.Cmd.onecmd(self, line)
        except AssertionError as err:
            if self.use_rawinput:
                print "\n".join("*** " + l for l in err.message.split('\n'))
            else:
                err.args = ('"{0}" => '.format(line) + err.message,)
                raise

    def do_exit(self, s):
        """Exit the interpreter. You can also use the Ctrl-D shortcut."""
        print "\n --- Exiting FlameTransfer. Bye!"
        sys.exit(0)
    do_EOF = do_exit
    do_quit = do_exit

    def read(self, prompt):
        if self.use_rawinput:
            try:
                line = raw_input(prompt)
            except EOFError:
                line = ''
        else:
            self.stdout.write(prompt)
            self.stdout.flush()
            line = self.stdin.readline()
            if not len(line):
                line = ''
            else:
                line = line.rstrip('\r\n')
        return line


class FlameTransferCmd(ShellCmd, SmartCmd, cmd.Cmd, object):
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
    flames = FlameDB()

    def precmd(self, line):
        """Reprint the line to know what is executed"""
        #print "\n >>> executing: ", line
        return line

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
                avbp_field     disc         cylinder         sphere 

                Awesome, right?
                """)

    def emptyline(self):
        """Empty line behavior: do nothing"""
        pass

    generators = "avbp_field disc cylinder rectangle brick sphere".split()
    def do_generate(self, s):
        args = s.split()
        assert len(args) == 2, "invalid number of arguments"
        typ, name = args
        typ = typ[0].lower()
        if typ == "d": # disc
            center = input_floats(self.read('Disk center (x y) : '))
            radius = input_floats(self.read('Disk radius (r)   : '))[0]
            self.flames.current = ActiveFlame(name, self.hip_exec)
            self.flames.current.define_flame_disc(center, radius)
        elif typ == "r": # rectangle
            pt_min = input_floats(self.read('Min point  (x y) : '))
            pt_max = input_floats(self.read('Max point  (x y) : '))
            self.flames.current = ActiveFlame(name, self.hip_exec)
            self.flames.current.define_flame_rectangle(pt_min, pt_max)
        elif typ == "b": # brick
            pt_min = input_floats(self.read('Min point  (x y z) : '))
            pt_max = input_floats(self.read('Max point  (x y z) : '))
            self.flames.current = ActiveFlame(name, self.hip_exec)
            self.flames.current.define_flame_brick(pt_min, pt_max)
        elif typ == "s": # sphere
            center = input_floats(self.read('Sphere center (x y z) : '))
            radius = input_floats(self.read('Sphere radius (r)     : '))
            self.flames.current = ActiveFlame(name, self.hip_exec)
            self.flames.current.define_flame_sphere(center, radius)
        elif typ == "c": # cylinder
            center = input_floats(self.read('Cylinder center (x y z) : '))
            radius = input_floats(self.read('Cylinder radius (r)     : '))
            vector = input_floats(self.read('Cylinder vector (x y z) : '))
            self.flames.current = ActiveFlame(name, self.hip_exec)
            self.flames.current.define_flame_cylinder(center, radius, vector)
        elif typ == "a": # avbp field
            avbp_mesh = self.read('AVBP mesh : ')
            assert isfile(avbp_mesh), "not such file"
            avbp_sol = self.read('AVBP sol : ')
            assert isfile(avbp_sol); "not such file"
            avbp_field = self.read('Scalar field : ')
            avbp_thresh = self.read('Threshold value : ')
            self.flames.current = ActiveFlame(name, self.hip_exec)
            self.flames.current.define_threshold_flame(
                    avbp_mesh, avbp_sol, avbp_field, avbp_thresh)
        else:
            raise AssertionError("unkown flame shape")
        self.flames.append(self.flames.current)

    def help_generate(self):
        print dedent("""\
                Generate an analytical flame shape. All lenths are in [m]
                > ge(nerate) <shape> <name> 
                |  <shape> : one of {}
                |  <name>  : careful, same name flame is overwritten""".format(
                  ", ".join(self.generators)))

    def complete_generate(self, text, line, begidx, endidx):
        return [ f for f in self.generators if f.startswith(text) ]

    writers = "mesh flame full".split()
    def do_write(self, s):
        assert len(self.flames) > 0, "no flames defined yet"
        assert self.flames.current.metas.pt_ref is not None, (
                "please set reference point and vector before writing")
        if s[:2] == "me": # me(sh)
            [flame.make_mesh() for flame in self.flames]
        elif s[:2] == "fl": # fl(ame)
            [flame.write_h5() for flame in self.flames]
        elif s[:2] == "fu": # fu(ll)
            [flame.write_full() for flame in self.flames]
        else:
            raise AssertionError("unknown option")

    def help_write(self):
        print dedent("""\
                Write output files
                > wr(ite) me(sh)|fl(ame)|nt(au)
                |  mesh  : write hdf5 mesh file
                |  flame : write hdf5 flame file
                |  full  : write hdf5 mesh, flame and matching xmf file""")

    def complete_write(self, text, line, begidx, endidx):
        return [ f for f in self.writers if f.startswith(text) ]

    listers = "flames refs metas json dump".split()
    def do_list(self, s):
        assert len(self.flames) > 0, "no flames defined yet"
        out = open(s.split()[-1], 'w') if (len(s.split()) > 1) else sys.stdout
        if out == sys.stdout: print
        arg = s[0].lower()
        if arg == "f": # flames
            width = max([len(f.metas.name) for f in self.flames])
            out.write("Cu   Name" + " "*width + "ndim  Volume     PtRef           VecRef\n")
            out.write("-" * (48 + width) + "\n")
            for i, f in enumerate(self.flames) :
                star = "*" if f == self.flames.current else " "
                out.write("{0} {1:2} {2:{width}}    {3}   {4:.2E}   {5: <15} {6}\n".format(
                    star, i+1, f.metas.name, f.metas.ndim, f.metas.volume,
                    f.metas.pt_ref, f.metas.vec_ref, width=width+2))
        elif arg == "r": # refs
            out.write("Reference point  : {}\n".format(self.flames.current.metas.pt_ref))
            out.write("Reference vector : {}\n".format(self.flames.current.metas.vec_ref))
        elif arg == "m": # metas
            out.write(" Key                  | Value\n")
            out.write(" -------------------------------------------\n")
            for k, v in sorted(self.flames.current.metas.static.items()):
                out.write(" {0:20} | {1}".format(k, v).replace('\n', ' ') + '\n')
            for k, v in sorted(self.flames.current.metas.vects.items()):
                out.write(" {0:20} | {1}".format(
                    k, v() if v is not None else None).replace('\n', ' ') + '\n')
        elif arg == "j": # json formatted output
            try:
                for flame in self.flames:
                    with open("{}.metas.json".format(flame.metas.name),
                              'w') as f:
                        flame.metas.write_json(f)
            except ImportError:
                raise AssertionError("could not import python module json")
        elif arg == "d": # dump
            assert out != sys.stdout, "<key> argument mandatory for dump"
            key = s.split()[-1]
            assert key in self.flames.current.metas.__dict__.keys(), (
                    "unknown key " + key)
            out.write("{}".format(
                self.flames.current.metas.__dict__[key]).translate(None, '[],'))
        else:
            raise AssertionError("unknown argument")
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

    setters = "current meta refs ntau static vector ".split()
    def do_set(self, s):
        """Change flame meta parameter"""
        assert len(self.flames) > 0, "no flames defined yet"
        args = s.split()
        switch = s[0].lower()
        if switch == "c": # current
            assert len(args) == 2, "invalid number of arguments"
            fla_nb, flame = self._flame_finder(args[1])
            self.flames.current = self.flames[fla_nb]
        elif switch in "sm": # static meta | mutable meta
            assert len(args) == 3, "invalid number of arguments"
            name, typ = args[1:]
            typ = typ[0].lower()
            if switch == "s":
                typ_dict = {'i': int,
                            'f': float,
                            'a': input_floats,
                            's': lambda s: s,} # there is no 'identity' func in python!
                assert typ in typ_dict.keys(), "unknown type " + typ
                val = typ_dict[typ](self.read("Value : "))
                setattr(self.flames.current.metas, name, val)
            else:
                typ_dict = {'n': NormalVector,
                            'v': Vector,
                            'p': Point,}
                assert typ in typ_dict.keys(), "unknown type " + typ
                val = typ_dict[typ](input_floats(self.read("Value : ")))
                self.flames.current.metas.set_vect(name, val)
        elif switch == "r": # refs
            assert len(args) == 1, "invalid number of arguments"
            pt_ref = input_floats(self.read('Reference point  (x y z) : '))
            vec_ref = input_floats(self.read('Reference vector  (x y z) : '))
            self.flames.current.metas.pt_ref = Point(pt_ref)
            self.flames.current.metas.vec_ref = NormalVector(vec_ref)
        elif switch == "n": # ntau
            assert len(args) == 1, "invalid number of arguments"
            typ = self.read("Type of N : ").lower()[:2]
            assert typ in "n1 n2 n3 cr gl ad".split(), "unknown N type"
            path = self.read("Path to N-tau file : ")
            assert isfile(path), "file can't be accessed"
            def exp_wrapper(path):
                with open(path, 'r') as f:
                    return StringIO(f.read().lower().replace('d', 'e'))
            freq, n, tau = np.loadtxt(exp_wrapper(path)).T
            typ = typ[:2].lower()
            if typ in "n1 cr".split():
                area, p_mean, gamma = input_floats(self.read("Area, Pmean, gamma : "))
                self.flames.current.set_n1_tau(freq, tau, n, area, p_mean, gamma)
            elif typ in "n2 gl".split():
                self.flames.current.set_n2_tau(freq, tau, n)
            elif typ in "n3 ad".split():
                u_bar, q_bar = input_floats(self.read("U_bar, Q_bar : "))
                self.flames.current.set_n3_tau(freq, tau, n, u_bar, q_bar)
        elif switch == "o": # order swap
            assert len(args) == 3, "invalid number of arguments"
            nb1, fl1 = self._flame_finder(args[1])
            nb2, fl2 = self._flame_finder(args[2])
            self.flames[nb1] = fl2
            self.flames[nb2] = fl1
        else:
            raise AssertionError("unknown settable")

    def help_set(self):
        print dedent("""\
                Set a value or parameter for the current flame
                > set current <ref>
                |  Set current flame to <ref>
                |  <ref> = either a flame name or number
                |  See ft> list flames to view which is current
                > set static|mutable meta_name type
                |  Set a meta parameter for the flame
                |  static types: string, int, float or array
                |  mutable types: normal_vector, vector, point
                |  Mutable metas will be affected by transformations (rotation...)
                > set refs
                |  Set the reference point and vector
                > set ntau
                |  Set values of N and tau
                |  Accepted N types: {n1|crocco, n2|global, n3|adim}
                |  See ?N_tau for more
                > set order fl1 fl2
                |  Swap 2 flames in list (by number or name)""")

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
        assert len(args) in [1, 2], "invalid number of arguments"
        assert isfile(args[0]), "no such file"
        flame = ActiveFlame('tmp', self.hip_exec)
        if args[-1][0] == "m":
            flame.load_metas(args[0])
        else:
            flame.load(args[0])
        self.flames.append(flame)
        self.flames.current = self.flames[-1]

    def help_read(self):
        print dedent("""\
                Read hdf5 flame file
                > re(ad) <path> [metas_only]
                |  <path>       : absolute or relative path to flame .h5 file
                |  [metas_only] : only read the flame's metas. Faster, but incomplete""")

    def complete_read(self, text, line, begidx, endidx):
        candidates = [f for f in os.listdir('.') if f[-9:] == ".flame.h5"]
        return [f for f in candidates if f.startswith(text)]

    def do_copy(self, s):
        assert len(self.flames) > 0, "no flames defined yet"
        args = s.split()
        assert len(args) == 2, "invalid number of arguments"
        src, dest = args
        assert not is_number(dest), "flame name cannot be a number"
        assert dest not in [flame.name for flame in self.flames], (
                "{} is already a declared flame".format(dest))
        fla_nb, flame = self._flame_finder(src)
        self.flames.append(copy.deepcopy(self.flames[fla_nb]))
        self.flames[-1].metas.name = dest
        self.flames.current = self.flames[-1]

    def help_copy(self):
        print dedent("""\
                Copy existing flame
                > co(py) <name_src> <name_dest>
                |  <name_src>  : name or number of flame to copy
                |  <name_dest> : name of destination flame""")

    transforms = "translate scale rotate".split()
    def do_transform(self, s):
        """Appy transformation to current flame"""
        assert len(self.flames) > 0, "no flames defined yet"
        switch = s[0].lower()
        if switch == "t":
            vect = input_floats(self.read('Translation vector : '))
            self.flames.current.transform(switch, vect)
        elif switch == "s":
            vect = input_floats(self.read('Scale vector : '))
            self.flames.current.transform(switch, vect)
        elif switch == "r":
            axis = self.read("Rotation axis : ")
            angle = self.read("Rotation angle (degrees) : ")
            self.flames.current.transform(switch, axis, float(angle))
        else:
            raise AssertionError("accepted transforms : " + " ".join(self.transforms))

    def help_transform(self):
        print dedent("""\
                Apply transformation to current flame
                > transform translate|scale|rotate""")

    def complete_transform(self, text, line, begidx, endidx):
        return [ f for f in self.transforms if f.startswith(text) ]

    exports = "AVSP5.6".split()
    def do_export(self, s):
        """Export current flame"""
        assert len(self.flames) > 0, "no flames defined yet"
        if s.lower() == "avsp5.6":
            avsp_mesh = self.read('AVSP mesh file : ')
            avsp_sol = self.read('AVSP sol file : ')
            self.flames.export(avsp_mesh, avsp_sol, "avsp5.6")
        else:
            raise AssertionError("unknown export format")

    def help_export(self):
        print dedent("""\
                Export all flames to another format
                > export AVSP5.6
                |  variable frequency n-tau - multi-flame AVSP version""")

    def complete_export(self, text, line, begidx, endidx):
        return [ f for f in self.exports if f.startswith(text) ]

    imports = "AVSP5.6".split()
    def do_import(self, s):
        """Import flames from other format"""
        if s.lower() == "avsp5.6":
            avsp_sol = self.read('AVSP sol file : ')
            self.flames.import_flames(avsp_sol, "avsp5.6", self.hip_exec)
        else:
            raise AssertionError("unknown import format")

    def help_import(self):
        print dedent("""\
                Import all flames from another format
                > import AVSP5.6
                |  variable frequency n-tau - multi-flame AVSP version""")

    def do_echo(self, s):
        """Print something on the command line. For batch debugging"""
        print '"' + s + '"'

    def complete_import(self, text, line, begidx, endidx):
        return [ f for f in self.imports if f.startswith(text) ]

    def _flame_finder(self, key):
        """Find a flame using it's name or it's number"""
        if is_number(key):
            nb = int(key)
            assert nb <= len(self.flames), (
                    "there are only {} flames declared".format(len(self.flames)))
            return nb-1, self.flames[nb-1]
        else:
            for i, fla in enumerate(self.flames):
                if fla.metas.name == key:
                    return i, fla
            raise AssertionError("no such flame declared")


if __name__ == '__main__':
    stdin = None if len(sys.argv) == 1 else open(sys.argv[-1])
    interpreter = FlameTransferCmd(stdin=stdin)
    if stdin is not None:
        interpreter.use_rawinput = False
    interpreter.cmdloop_with_keyboard_interrupt()
