# FlameTransfer

Welcome to the FlameTransfer app git repository

## Installation instructions

FlameTransfer is a python tool.  As such, no compilation is necessary: simply
clone the git repository somewhere appropriate, and execute it.  It does however
rely on a short list of dependancies:

 - _hip:_ FlameTransfer is in fact largely an encapsulation of many hip
   commands that could be performed manually.  The environment variable
   `$HIPEXEC` must be set to an up to date hip executable, *e.g*. at CERFACS
   on the NFS:

   ```bash
   export HIPEXEC="/home/cfd2/avbp/HIP/latest/hip-latest"
   ```

 - **python 2.7 and some modules**: _numpy_, _h5py_.  The usual installations
   using pip are perfect for this job, if you have root access.  Otherwise,
   find a machine that has these installed.  To check this, you can try *e.g.*:

   ```bash
   python -c "import numpy; import h5py; print 'All good'"
   ```

## Usage

To run, simply type

```bash
python /path/to/flametransfer.py
```

or add it to your `$PATH`.  Hopefully, the builtin help system is all you need
from here.

## Authors, Contribution

FlameTransfer tool was initially written by [Corentin J.
Lapeyre](mailto:lapeyre@cerfacs.fr).   It relies heavily on hip by Jens-Dominik
Mueller and Gabriel Staffelbach.  Please send any feedback to Corentin.

