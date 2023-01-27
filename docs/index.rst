.. _home:

.. figure:: /_static/images/NuMAD-header.png
   :target: https://github.com/sandialabs/NuMAD


##################################################################################
Numerical Manufacturing And Design (NuMAD) Tool for Wind Turbine Blades
##################################################################################

.. toctree::
   :maxdepth: 2
   :hidden:
   
   Home <self>
   API documentation <apidoc/pynumad>
   introduction/index.rst
   user/index.rst
   developer.rst
   gettingstarted





The structural design and optimization of wind turbine blades is a
complex task. In many cases it is difficult to find the optimal design
of a turbine blade by hand, or by trial and error, and the software
tools used for such designs are most effective when integrated into
automated optimization and analysis algorithms. A new version of the
software tool `NuMAD (Numerical Manufacturing And Design) <https://github.com/sandialabs/NuMAD>`_ for the design
and modeling of wind turbine blades is developed and described. 

Newly released :ref:`NuMADv3` is structured to be run from a scripting environment
and easily called by optimization processes. 
The previous release, :ref:`NuMADv2`, relied on the use of a graphical
user interface. Modifications were made to the MATLAB-based source code
to decouple the internal functions from the graphical user interface,
and additional functionality was added for convenience in performing
high-fidelity finite element analysis. NuMAD 3.0 has been successfully
implemented for optimization of large, flexible rotor blades, and is
suitable to be made available for public use.


.. _developers:

NuMAD Developers
=====================
NuMAD has been developed by `Sandia National Laboratories 
(Sandia) <https://energy.sandia.gov/programs/renewable-energy/wind-power/>`_, 
funded by the U.S. Department of Energy’s Wind Energy Technologies Technologies Office. 


Current members of the development team include:

* Joshua Paquette (“Josh”) (Sandia - PI)
* Evan Anderson (Sandia)
* Ernesto Camarena (Sandia)
* Ryan James Clarke (Sandia)
* Kirk Bonney (Sandia)

Former members of the development team include:

* Jonathan Charles Berg
* Brian Resor
* Daniel Laird
* Kelley Ruehl (Sandia)
* Christopher Lee Kelley (Sandia)
* Brandon Lee Ennis (Sandia)

Funding
=======

Development and maintenance of the NuMAD code is funded by the U.S. Department of Energy’s Wind Energy Technologies Office.

Sandia National Laboratories is a multi-mission laboratory managed and operated by National Technology and Engineering Solutions of Sandia, LLC., a wholly owned subsidiary of Honeywell International, Inc., for the U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
