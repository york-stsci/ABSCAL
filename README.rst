Absolute Flux calibration software converted from Dr. Ralph Bohlin's IDL code.
------------------------------------------------------------------------------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

STIPS and WFC3 calibration software preserving the functionality of Ralph
Bohlin's IDL code in python.

*Ensuring that 30 years of calibration software are useful for the next 30 years.*

*Scio enim passus sum, et nunc vicissim tuum est*

Installation and Documentation
------------------------------

Because the ABSCAL documentation is currently only available online in source form, a 
quick cheat sheet on installing ABSCAL and compiling its documentation is provided below:

1. Clone the ABSCAL source from the `github repository <https://github.com/spacetelescope/ABSCAL>`_::

    git clone https://github.com/spacetelescope/ABSCAL.git

    cd ABSCAL

2. Use the `env_setup.py` script to set up the ABSCAL environment::

    python env_setup.py

    conda activate abscal

4. Build the HTML documentation (this has not been tested recently)::

    cd docs
    
    make html

The documentation webpages will now be in `ABSCAL/docs/_build/html`.

License
-------

This project is Copyright (c) Space Telescope Science Institute and licensed under
the terms of the BSD 3-Clause license.


Contributing
------------

We love contributions! abscal is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
abscal based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.
