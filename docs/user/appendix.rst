.. _appendix:


Appendix
========

.. _bladeClass:

Blade Class
------------------

.. autoclass:: numadObjects.BladeDef
	:members: 
	:exclude-members: 
	:no-undoc-members: 
	
.. TODO: properties and methods should be autopopulated (above), manual for now (below)

Blade Properties	
~~~~~~~~~~~~~~~~~~

.. autoattribute:: numadObjects.BladeDef.aerocenter
.. autoattribute:: numadObjects.BladeDef.chord
.. autoattribute:: numadObjects.BladeDef.chordoffset
.. autoattribute:: numadObjects.BladeDef.components
.. autoattribute:: numadObjects.BladeDef.degreestwist
.. autoattribute:: numadObjects.BladeDef.ispan
.. autoattribute:: numadObjects.BladeDef.leband
.. autoattribute:: numadObjects.BladeDef.materials
.. autoattribute:: numadObjects.BladeDef.mesh
.. autoattribute:: numadObjects.BladeDef.naturaloffset
.. autoattribute:: numadObjects.BladeDef.percentthick
.. autoattribute:: numadObjects.BladeDef.prebend
.. autoattribute:: numadObjects.BladeDef.rotorspin
.. autoattribute:: numadObjects.BladeDef.span
.. autoattribute:: numadObjects.BladeDef.sparcapoffset
.. autoattribute:: numadObjects.BladeDef.sparcapwidth
.. autoattribute:: numadObjects.BladeDef.stations
.. autoattribute:: numadObjects.BladeDef.sweep
.. autoattribute:: numadObjects.BladeDef.swtwisted
.. autoattribute:: numadObjects.BladeDef.teband


Blade Methods
~~~~~~~~~~~~~~~~~~
.. NOTE:`set.naturaloffset` syntax isn't compatible with autodocs, e.g. `setNaturaloffset` would be
.. TODO: .. autoattribute:: numadObjects.BladeDef.set.naturaloffset
.. TODO: .. autoattribute:: numadObjects.BladeDef.set.rotorspin
.. TODO: .. autoattribute:: numadObjects.BladeDef.set.swtwisted
.. TODO: .. autoattribute:: numadObjects.BladeDef.findLayerExtents

.. autoattribute:: numadObjects.BladeDef.addStation
.. autoattribute:: numadObjects.BladeDef.addComponent
.. autoattribute:: numadObjects.BladeDef.addMaterial
.. autoattribute:: numadObjects.BladeDef.updateBlade
.. autoattribute:: numadObjects.BladeDef.updateGeometry
.. autoattribute:: numadObjects.BladeDef.updateKeypoints
.. autoattribute:: numadObjects.BladeDef.updateBOM
.. autoattribute:: numadObjects.BladeDef.readYAML
.. autoattribute:: numadObjects.BladeDef.writeYAML
.. autoattribute:: numadObjects.BladeDef.generateBeamModel
.. autoattribute:: numadObjects.BladeDef.generateFEA 
.. autoattribute:: numadObjects.BladeDef.writeBOMxls
.. autoattribute:: numadObjects.BladeDef.writePlot3D
.. autoattribute:: numadObjects.BladeDef.getprofileTEtype
.. autoattribute:: numadObjects.BladeDef.downsampleProfile
.. autoattribute:: numadObjects.BladeDef.delete
.. autoattribute:: numadObjects.BladeDef.surf
.. autoattribute:: numadObjects.BladeDef.plotregions
.. autoattribute:: numadObjects.BladeDef.plotgeom
.. autoattribute:: numadObjects.BladeDef.plotbom
.. autoattribute:: numadObjects.BladeDef.plotprofile

.. NOTE: these functions are not included within the object, they are appended at the end
.. TODO: .. autoattribute:: numadObjects.BladeDef.findLayerExtents
.. TODO: .. autoattribute:: numadObjects.BladeDef.findRegionExtents
.. TODO: .. autoattribute:: numadObjects.BladeDef.getTEtype
.. TODO: .. autoattribute:: numadObjects.BladeDef.fprintf_matrix



.. Kelley: remove these tables
.. .. _bladeDefTable:
.. .. csv-table:: ``BladeDef``: A class definition for wind & water turbine blades.
..    :file: bladeDefTable.csv
..    :widths: 1,2
..    :header-rows: 1


.. .. _bladeDefMethodsTable:
.. .. csv-table::  List of Methods for ``BladeDef`` class.
..    :file: bladeDefMethodsTable.csv
..    :widths: 1
..    :header-rows: 1


.. _materialClass:

Material Class
------------------

.. autoclass:: numadObjects.MaterialDef
	:members: 
	:exclude-members: 
	:no-undoc-members: 

.. TODO: properties and methods should be autopopulated (above), manual for now (below)	
	
Material Properties	
~~~~~~~~~~~~~~~~~~~

.. autoattribute:: numadObjects.MaterialDef.name
.. autoattribute:: numadObjects.MaterialDef.type
.. autoattribute:: numadObjects.MaterialDef.layerthickness
.. autoattribute:: numadObjects.MaterialDef.ex
.. autoattribute:: numadObjects.MaterialDef.ey
.. autoattribute:: numadObjects.MaterialDef.ez
.. autoattribute:: numadObjects.MaterialDef.gxy
.. autoattribute:: numadObjects.MaterialDef.gyz
.. autoattribute:: numadObjects.MaterialDef.gxz
.. autoattribute:: numadObjects.MaterialDef.prxy
.. autoattribute:: numadObjects.MaterialDef.pryz
.. autoattribute:: numadObjects.MaterialDef.prxz
.. autoattribute:: numadObjects.MaterialDef.density
.. autoattribute:: numadObjects.MaterialDef.drydensity
.. autoattribute:: numadObjects.MaterialDef.uts
.. autoattribute:: numadObjects.MaterialDef.ucs
.. autoattribute:: numadObjects.MaterialDef.uss
.. autoattribute:: numadObjects.MaterialDef.xzit
.. autoattribute:: numadObjects.MaterialDef.xzic
.. autoattribute:: numadObjects.MaterialDef.yzit
.. autoattribute:: numadObjects.MaterialDef.yzic
.. autoattribute:: numadObjects.MaterialDef.g1g2
.. autoattribute:: numadObjects.MaterialDef.alp0
.. autoattribute:: numadObjects.MaterialDef.etat
.. autoattribute:: numadObjects.MaterialDef.etal
.. autoattribute:: numadObjects.MaterialDef.m
.. autoattribute:: numadObjects.MaterialDef.reference
	
.. Kelley: remove these tables
   
.. _materialDefTable:
.. .. csv-table:: ``MaterialDef``: A class definition for blade materials. Materials properties are defined in the principal material coordinate system.
..    :file: materialDefTable.csv
..    :widths: 1,2
..    :header-rows: 1


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


   
   

