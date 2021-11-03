.. _gettingStarted:

Getting Started 
================

Software Requirements
----------------------
NuMAD is developed in MATLAB, and requires the following software packages:


==========================  						=============================
**Required Packages**        						**Oldest Compatible Version**
`MATLAB <https://www.mathworks.com/products/matlab.html>`_  		Version 9.9  (R2020b)
`FAST (Optional) <https://www.nrel.gov/wind/nwtc/fastv7.html>`_ 	v7.02.00d
`ANSYS (Optional) <https://www.ansys.com/>`_	    			R12
`PreComp (Optional) <https://www.nrel.gov/wind/nwtc/precomp.html>`_ 	v1.00.03
`TurbSim (Optional) <https://www.nrel.gov/wind/nwtc/turbsim.html>`_     v1.50
`IECWind (Optional) <https://www.nrel.gov/wind/nwtc/iecwind.html>`_     v5.01.01
`BModes (Optional) <https://www.nrel.gov/wind/nwtc/bmodes.html>`_       v3.00.00
`Crunch (Optional) <https://www.nrel.gov/wind/nwtc/crunch.html>`_       v3.00.00
========================== 						=============================
.. Unsupported IECWind, ADAMS
 

Installing NuMAD 
----------------
The code containing all open-source tools for NuMAD can be downloaded or
cloned from the `NuMAD GitHub repository <https://github.com/sandialabs/NuMAD>`_. 

There are 3 directories within ``NuMAD/source`` : ``toolbox``,
``objects``, and ``optimization``, their purpose is described in the table below:  

============================ ===================================================
Source Directory       	 	Description
============================ ===================================================
``numadGUI``			contains files related to the graphical user interface (GUI) which was introduced in NuMAD v2.0.
``objects``			contains files associated to the class definition of the blade object and ``.yaml`` file capabilities.
``optimization``		contains bundles of tools for several purposes, including ``runIEC``, as explained further in the :ref:`AeroelasticSimRunIEC` section, file processing functions for input and output from other programs such as FAST and Crunch, and setup and execution of ANSYS models for analyzing quantities like material rupture, fatigue, and buckling under loading
============================ ===================================================


Installation Steps
~~~~~~~~~~~~~~~~~~

1.    Install `MATLAB <https://www.mathworks.com/products/matlab.html>`_
2.    Open ``addNumadPaths.m`` in an editor and set the ``numadPath`` to point to the source directory in the repo.
3.    If you plan to use ANSYS or other analysis programs, update the corresponding path variable to the actual path of each executable you plan to use. Save and close ``addNumadPaths.m.``
4.    Move the ``addNumadPaths.m`` file to the MATLAB directory of your user account (e.g. ``C:\Users\username\Documents\MATLAB``)

Every time a new MATLAB session is started type ``addNumadPaths`` to be able to run NuMAD.

