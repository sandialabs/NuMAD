.. _gettingStarted:

Getting Started 
================

Software Requirements
----------------------
NuMAD is developed in MATLAB, and requires the following software packages:

==========================  =============================
**Required Packages**        **Oldest Compatible Version**
MATLAB                      Version 9.9  (R2020b)
FAST (Optional)		    v7
ANSYS (Optional)	    R12
PreComp (Optional)      v1.00.03
BModes (Optional)       v3.00.00
Crunch (Optional)       v3.00.00
==========================  =============================

 

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


Installation Steps
~~~~~~~~~~~~~~~~~~

1.    Install `MATLAB <https://www.mathworks.com/products/matlab.html>`_
2.    Open ``addNumadPaths.m`` in an editor and set the ``numadPath`` to point to the source directory in the repo.
3.    If you plan to use ANSYS or other analysis programs, update the corresponding path variable to the actual path of each executable you plan to use. Save and close ``addNumadPaths.m.``
4.    Move the ``addNumadPaths.m`` file to the MATLAB directory of your user account (e.g. ``C:\Users\username\Documents\MATLAB``)

Every time a new MATLAB session is started type ``addNumadPaths`` to be able to run NuMAD.

