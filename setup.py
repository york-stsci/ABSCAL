from setuptools import setup
from setuptools import find_packages

from glob import glob, iglob
import os


file_dir = os.path.abspath(__file__)
version_str = '2.0.dev'

setup(
    name = 'abscal',
    description = 'HST WFC3 and STIS absolute flux calibration',
    url = 'https://github.com/spacetelescope/ABSCAL',
    author = 'Brian York, Ralph Bohlin, Susana Deustua, Karl Gordon',
    author_email = 'york@stsci.edu, bohlin@stsci.edu, deustua@stsci.edu, kgordon@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python'],
    packages = find_packages(),
    install_requires = [
                        "numpy", 
                        "scipy", 
                        "astropy>=3", 
                        "photutils",
                        "matplotlib",
                        "rebin",
                        "stistools",
                        # Removed because it's not on pip
#                         "hstcal",
                        "astroquery"
                       ],
    version = version_str,
    scripts=glob("abscal/commands/*"),
    package_data =  {
                        "": [x.replace("abscal/", "") for x in iglob('abscal/data/**/*', recursive=True)] +
                            [x.replace("abscal/", "") for x in iglob('abscal/idl_code/**/*', recursive=True)],
#                             ["data/*", 
#                              "data/defaults/*", 
#                              "idl_code", 
#                              "idl_code/*",
#                              "idl_code/common/*", 
#                              "idl_code/cookbooks/*", 
#                              "idl_code/stis/*",
#                              "idl_code/stis/calstis/*", 
#                              "idl_code/stis/manual/*",
#                              "idl_code/wfc3/*"],
                        "wfc3": [x.replace("abscal/wfc3/", "") for x in iglob('abscal/wfc3/data/**/*', recursive=True)]
#                                  "data/*",
#                                  "data/pnref/*",
#                                  "data/wfcref/*"]
                    },
    include_package_data=True
    )
