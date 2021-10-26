.. _gettingStarted:

Getting Started 
================

Software Requirements
----------------------
NuMAD is developed in MATLAB, and requires the following software packages:

==========================  =============================
**Required Toolbox**        **Oldest Compatible Version**
MATLAB                      Version 9.9  (R2020b)
FAST (Optional)		    v7
ANSYS (Optional)
PreComp (Optional)
BModes (Optional)
FAST (Optional)
PLOT3D (Optional)
==========================  =============================

.. Kelley: add version requirements
 

Installing NuMAD 
----------------
The code containing all open-source tools for NuMAD can be downloaded or
cloned from the `NuMAD GitHub repository <https://github.com/sandialabs/NuMAD>`_. 
There are 3 directories within ``NuMAD/source`` : ``toolbox``,
``objects``, and ``optimization``, their purpose is described in the table below:  

============================ ===================================================
Source Directory       	 	Description
============================ ===================================================
``toolbox``			contains basic functions and operations needed for performing analysis with packages such as precomp, BPE, and ANSYS
``objects``			contains the class definition of the blade object, which stores the geometric, airfoil and material data for a given blade design
``optimization``		contains bundles of tools for several purposes, including ``runIEC``, as explained further in the :ref:`AeroelasticSimRunIEC` section, file processing functions for input and output from other programs such as FAST and Crunch, and setup and execution of ANSYS models for analyzing quantities like material rupture, fatigue, and buckling under loading
============================ ===================================================

It is advisable to save the NuMAD source
in a directory called ``DesignCodes\NuMAD``, for example ``C:\DesignCodes\NuMAD\source\``. 
Before beginning any design or analysis of a
blade, the paths of these folders should be added to the MATLAB domain
of working directories, using the script ``addNumadPaths.m``. All tools and
functions distrubuted with the NuMAD source code will then be accessible from from a MATLAB
script or the command window. If you intend to use NuMAD in conjuction with ANSYS then the global variable 
``ANSYS_Path`` needs to be modified to your ANSYS working directory file path name.

Installation Steps
~~~~~~~~~~~~~~~~~~

1.    Install `MATLAB <https://www.mathworks.com/products/matlab.html>`_
2.    Open ``addNumadPaths.m`` in an editor and set the ``NuMAD_path`` to the source directory in the repo. It is highly recommended that the contents of the repo be located at ``C:\DesignCodes\NuMAD``.
3.    If you plan to use ANSYS, update the ``ANSYS_Path`` variable to the actual path to your ANSYS executable. Save and close ``addNumadPaths.m.``
4.    Move the ``addNumadPaths.m`` file to the MATLAB directory of your user account (e.g. ``C:\Users\username\Documents\MATLAB``)


.. Note::
	If a path other than ``C:\DesignCodes\NuMAD\source\`` is used, path definitions in the ``runIEC_ipt.m`` input settings script may need to be modified, as described further in the :ref:`AeroelasticSimRunIEC` section, and the :ref:`appendix`. 