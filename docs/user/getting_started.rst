.. _gettingStarted:

Getting Started 
================

MATLAB Requirements
-------------------

NuMAD is developed in MATLAB, and requires the following toolboxes:

.. Kelley: add MATLAB version requirements
 
==========================  =============================
**Required Toolbox**        **Oldest Compatible Version**
MATLAB                      Version 9.9  (R2020b)
==========================  =============================


Installing NuMAD 
----------------
The code containing all open-source tools for NuMAD can be downloaded or
cloned from the `NuMAD GitHub repository <https://github.com/sandialabs/NuMAD>`__. 
There are 3 main folders within the ``NuMAD/source`` directory: ``NuMAD_toolbox``,
``preNuMAD``, and ``rotor_optimization``. It is advisable to save the NumAD source
in a directory called ``DesignCodes`` on, for example ``C:\DesignCodes\NuMAD\source\``. 


Before beginning any design or analysis of a
blade, the paths of these folders should be added to the MATLAB domain
of working directories, using the script ``paths.m``. All tools and
functions distrubuted with the NuMAD source code will then be accessible from from a MATLAB
script or the command window. If you intend to use NuMAD in conjuction with ANSYS then the global variable 
``ANSYS_Path`` needs to be modified to your ANSYS working directory file path name.

============================ ===================================================
Source Directory       	 	Description
============================ ===================================================
``NuMAD_toolbox``		contains basic functions and operations needed for performing analysis with packages such as precomp, BPE, and ANSYS
``preNuMAD``			contains the class definition of the blade object, which stores the geometric, airfoil and material data for a given blade design
``rotor_optimization``		contains bundles of tools for several purposes, including ``runIEC``, as explained further in the :ref:`AeroelasticSimRunIEC` section, file processing functions for input and output from other programs such as FAST and Crunch, and setup and execution of ANSYS models for analyzing quantities like material rupture, fatigue, and buckling under loading
============================ ===================================================




.. Note::
	If a path other than ``C:\DesignCodes\NuMAD\source\`` is used, path definitions in the ``runIEC_ipt.m`` input settings script may need to be modified, as described further in the :ref:`AeroelasticSimRunIEC` section, and the :ref:`appendix`. 