.. _appendix:


Appendix
========

Blade Class
------------------

.. _bladeDefTable:
.. csv-table:: ``BladeDef``: A class definition for wind & water turbine blades.
   :file: bladeDefTable.csv
   :widths: 1,2
   :header-rows: 1


.. _bladeDefMethodsTable:
.. csv-table::  List of Methods for ``BladeDef`` class.
   :file: bladeDefMethodsTable.csv
   :widths: 1
   :header-rows: 1


Material Class
------------------
   
.. _materialDefTable:
.. csv-table:: ``MaterialDef``: A class definition for blade materials. Materials properties are defined in the principal material coordinate system.
   :file: materialDefTable.csv
   :widths: 1,2
   :header-rows: 1


Station Class
------------------

.. _stationDefTable:
.. csv-table::  ``StationDef``: A class definition for blade stations.
   :file: stationDefTable.csv
   :widths: 1,2
   :header-rows: 1


Component Class
------------------

.. _componentDefTable:
.. csv-table::  ``ComponentDef``: A class definition for blade components.
   :file: componentDefTable.csv
   :widths: 1,2
   :header-rows: 1


Airfoil Class
------------------

.. _airfoilDefTable:
.. csv-table::  ``AirfoilDef``: A class definition for airfoil profiles.
   :file: airfoilDefTable.csv
   :widths: 1,2
   :header-rows: 1
   

Stack Class
------------------ 

.. _stackDefTable:
.. table:: ``StackDef``: A class definition for a stack of composite layers.

    +---------------------+----------------------------------------------------+
    | **StackDef          | **Property Description**                           |
    | Property**          |                                                    |
    +=====================+====================================================+
    | ``name``            | Name of the stack / composite material used by     |
    |                     | NuMAD                                              |
    +---------------------+----------------------------------------------------+
    | ``plygroups``       | Array of ply structures, one for each ply.         |
    |                     |                                                    |
    |                     | ``ply = struct('component'),...% parent comp``     |
    |                     |                                                    |
    |                     | ``'materialid',[],...% materialid of ply``         |
    |                     |                                                    |
    |                     | ``'thickness',[],... % thickness [mm]``            |
    |                     |                                                    |
    |                     | ``% of single ply``                                |
    |                     |                                                    |
    |                     | ``'angle',[],... % ply angle``                     |
    |                     |                                                    |
    |                     | ``'nPlies',[]); % number of plies``                |
    +---------------------+----------------------------------------------------+
    | ``indices``         | Indices of stack                                   |
    |                     |                                                    |
    |                     | [in board station, out board station, 1\ :sup:`st` |
    |                     | kepoint, 2\ :sup:`nd` keypoint]                    |
    +---------------------+----------------------------------------------------+


   
   

