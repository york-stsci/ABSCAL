from setuptools import setup
from setuptools import find_packages

import glob
import os


file_dir = os.path.abspath(__file__)
version_str = '1.0.dev'

setup(
    name = 'abscal',
    description = 'HST WFC3 and STIS absolute flux calibration',
    url = 'https://github.com/spacetelescope/ABSCAL',
    author = 'Brian York, Ralph Bohlin, Susana Deustua',
    author_email = 'york@stsci.edu, bohlin@stsci.edu, deustua@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python'],
    packages = find_packages(),
    install_requires = [
                        "numpy", 
                        "scipy", 
                        "astropy>=3", 
                        "photutils",
                        "matplotlib",
                        "rebin"
                       ],
    version = version_str,
    scripts=glob.glob("abscal/commands/*"),
    package_data =  {
                        "": ["data/*", 
                             "data/defaults/*", 
                             "idl_code", 
                             "idl_code/*",
                             "idl_code/common/*", 
                             "idl_code/cookbooks/*", 
                             "idl_code/stis/*",
                             "idl_code/stis/calstis/*", 
                             "idl_code/stis/manual/*",
                             "idl_code/wfc3/*"],
                        "wfc3": ["data/*",
                                 "data/pnref/*",
                                 "data/wfcref/*"]
                    },
    include_package_data=True
    )
