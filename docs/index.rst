PyARAGORN |Stars|
=================

.. .. |Logo| image:: /_images/logo.png
..    :scale: 40%
..    :class: dark-light

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyaragorn.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyaragorn/stargazers
   :class: dark-light

*Cython bindings and Python interface to* `ARAGORN <https://www.trna.se/>`_, *a tRNA, tmRNA and mtRNA gene finder*.

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyaragorn/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyaragorn/actions
   :class: dark-light

.. |GitLabCI| image:: https://img.shields.io/gitlab/pipeline/larralde/pyaragorn/main?gitlab_url=https%3A%2F%2Fgit.embl.de&logo=gitlab&style=flat-square&maxAge=600
   :target: https://git.embl.de/larralde/pyaragorn/-/pipelines
   :class: dark-light

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyaragorn?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyaragorn/
   :class: dark-light

.. |PyPI| image:: https://img.shields.io/pypi/v/pyaragorn.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyaragorn
   :class: dark-light

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyaragorn?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyaragorn
   :class: dark-light

.. |AUR| image:: https://img.shields.io/aur/version/python-pyaragorn?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyaragorn
   :class: dark-light

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyaragorn?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyaragorn/#files
   :class: dark-light

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyaragorn.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyaragorn/#files
   :class: dark-light

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyaragorn.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyaragorn/#files
   :class: dark-light

.. |License| image:: https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/gpl-3.0/
   :class: dark-light

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/pyaragorn/
   :class: dark-light

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=3600
   :target: https://git.embl.de/larralde/pyaragorn/
   :class: dark-light

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyaragorn.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyaragorn/issues
   :class: dark-light

.. |Docs| image:: https://img.shields.io/readthedocs/pyaragorn?style=flat-square&maxAge=3600
   :target: http://pyaragorn.readthedocs.io/en/stable/?badge=stable
   :class: dark-light

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/pyaragorn/blob/main/CHANGELOG.md
   :class: dark-light

.. |Downloads| image:: https://img.shields.io/pypi/dm/pyaragorn?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/pyaragorn
   :class: dark-light


Overview
--------

PyARAGORN is a Python module that provides bindings to ARAGORN using
`Cython <https://cython.org/>`_. It directly interacts with the ARAGORN
internals, which has the following advantages:


.. grid:: 1 2 3 3
   :gutter: 1

   .. grid-item-card:: :fas:`battery-full` Batteries-included

      Just add ``pyaragorn`` as a ``pip`` dependency, no need
      for the ARAGORN binary.

   .. grid-item-card:: :fas:`screwdriver-wrench` Flexible I/O

      Directly pass sequences to process as Python `str` objects, no 
      need for intermediate files.

   .. grid-item-card:: :fas:`gears` Practical output

      Retrieve the results as `~pyaragorn.Gene` objects directly 
      without parsing output files.

   .. grid-item-card:: :fas:`check` Consistent results 

      Get the exact same results as ARAGORN ``1.2.41``.



Features
--------

This library wraps the original source code of ARAGORN ``1.2.41``, processed
with the `pycparser <https://github.com/eliben/pycparser>`_ library to perform 
AST transformation of the original code in order to eliminate global variables 
among other fixes. 


Setup
-----

Run ``pip install pyaragorn`` in a shell to download the latest release and 
its dependencies from PyPi, or have a look at the
:doc:`Installation page <guide/install>` to find other ways to 
install ``pyaragorn``.


Library
-------

Check the following pages of the user guide or the API reference for more
in-depth reference about library setup, usage, and rationale:

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. include:: related.rst


License
-------

This library is provided under the `GNU General Public License v3.0 <https://choosealicense.com/licenses/gpl-3.0/>`_
or later. ARAGORN and ARWEN were developed by Dean Laslett and are distributed 
under the terms of the GPLv3 or later as well. See the 
:doc:`Copyright Notice <guide/copyright>` section for the full GPLv3 license.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* `ARAGORN`_ *authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `Leiden University Medical Center <https://www.lumc.nl/>`_
*in the* `Zeller Lab <https://zellerlab.org>`_.

