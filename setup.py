#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name="FlameTransfer",
    version="1.0.0-rc2",
    packages=find_packages(),
    # scripts=['flametransfer/flametransfer.py'],
    entry_points={
        'console_scripts': [
            'flametransfer = flametransfer.flametransfer:main',
            'visu = flametransfer.visu:main'
        ],
    },

    install_requires=[
        'numpy>=1.10',
        'h5py>=2.5'
    ],

    # metadata
    author="Corentin J. Lapeyre",
    author_email="lapeyre@cerfacs.fr",
    description="Active flame manipulation tool",
    license="MIT",
    url="https://github.com/clapeyre/flametransfer",
)
