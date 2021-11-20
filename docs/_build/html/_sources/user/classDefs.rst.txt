.. _classDefs:

Object Classes, Properties, and Methods
=======================================

.. _bladeClass:

Blade Class
------------------
.. TODO: properties and methods should be autopopulated, manual for now (below)

.. autoclass:: numadObjects.BladeDef
	:members: 
	:exclude-members: 
	:no-undoc-members: 
	

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
.. autoattribute:: numadObjects.BladeDef.checkNaturalOffset
.. autoattribute:: numadObjects.BladeDef.checkRotorSpin
.. autoattribute:: numadObjects.BladeDef.checkSwtwisted
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
.. autoattribute:: numadObjects.BladeDef.findLayerExtents
.. autoattribute:: numadObjects.BladeDef.findRegionExtents
.. autoattribute:: numadObjects.BladeDef.getTEtype
.. autoattribute:: numadObjects.BladeDef.fprintf_matrix


.. _materialClass:

Material Class
------------------
.. TODO: properties and methods should be autopopulated, manual for now (below)	

.. autoclass:: numadObjects.MaterialDef
	:members: 
	:exclude-members: 
	:no-undoc-members: 

	
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


.. _StationClass:

Station Class
------------------
.. TODO: properties and methods should be autopopulated, manual for now (below)	

.. autoclass:: numadObjects.StationDef
	:members: 
	:exclude-members: 
	:no-undoc-members: 

	
Station Properties	
~~~~~~~~~~~~~~~~~~~
.. autoattribute:: numadObjects.StationDef.airfoil
.. autoattribute:: numadObjects.StationDef.spanlocation


.. _ComponentClass:

Component Class
------------------
.. TODO: properties and methods should be autopopulated, manual for now (below)	

.. autoclass:: numadObjects.ComponentDef
	:members: 
	:exclude-members: 
	:no-undoc-members: 

	
Component Properties	
~~~~~~~~~~~~~~~~~~~~
.. autoattribute:: numadObjects.ComponentDef.group
.. autoattribute:: numadObjects.ComponentDef.name
.. autoattribute:: numadObjects.ComponentDef.materialid
.. autoattribute:: numadObjects.ComponentDef.fabricangle
.. autoattribute:: numadObjects.ComponentDef.hpextents
.. autoattribute:: numadObjects.ComponentDef.lpextents
.. autoattribute:: numadObjects.ComponentDef.cp
.. autoattribute:: numadObjects.ComponentDef.imethod
.. autoattribute:: numadObjects.ComponentDef.pinnedends


.. _AirfoilClass:

Airfoil Class
------------------
.. TODO: properties and methods should be autopopulated, manual for now (below)	

.. autoclass:: numadObjects.AirfoilDef
	:members: 
	:exclude-members: 
	:no-undoc-members: 

	
Airfoil Properties	
~~~~~~~~~~~~~~~~~~~~
.. autoattribute:: numadObjects.AirfoilDef.name
.. autoattribute:: numadObjects.AirfoilDef.reference
.. autoattribute:: numadObjects.AirfoilDef.coordinates
.. autoattribute:: numadObjects.AirfoilDef.c
.. autoattribute:: numadObjects.AirfoilDef.camber	
.. autoattribute:: numadObjects.AirfoilDef.thickness
.. autoattribute:: numadObjects.AirfoilDef.percentthick
.. autoattribute:: numadObjects.AirfoilDef.maxthick
.. autoattribute:: numadObjects.AirfoilDef.TEtype
.. autoattribute:: numadObjects.AirfoilDef.x
.. autoattribute:: numadObjects.AirfoilDef.y


.. _StackClass:

Stack Class
------------------ 
.. TODO: properties and methods should be autopopulated, manual for now (below)	

.. autoclass:: numadObjects.StackDef
	:members: 
	:exclude-members: 
	:no-undoc-members: 

	
Stack Properties	
~~~~~~~~~~~~~~~~~~~~
.. autoattribute:: numadObjects.StackDef.name
.. autoattribute:: numadObjects.StackDef.plygroups
.. autoattribute:: numadObjects.StackDef.indices

Stack Methods
~~~~~~~~~~~~~~~~~~~~
.. autoattribute:: numadObjects.StackDef.addply



.. _IecClass:

IEC Class
------------------
.. TODO: properties and methods should be autopopulated, manual for now (below)

.. autoclass:: numadObjects.IECDef
	:members: 
	:exclude-members: 
	:no-undoc-members: 
	

IEC Properties	
~~~~~~~~~~~~~~~~~~
.. autoattribute:: numadObjects.IECDef.BldGagNd
.. autoattribute:: numadObjects.IECDef.Class
.. autoattribute:: numadObjects.IECDef.delay
.. autoattribute:: numadObjects.IECDef.designLife
.. autoattribute:: numadObjects.IECDef.fastsim
.. autoattribute:: numadObjects.IECDef.fatigueCriterion
.. autoattribute:: numadObjects.IECDef.fstfn
.. autoattribute:: numadObjects.IECDef.fullLoads
.. autoattribute:: numadObjects.IECDef.gageSetCase
.. autoattribute:: numadObjects.IECDef.lin
.. autoattribute:: numadObjects.IECDef.momentMaxRotation
.. autoattribute:: numadObjects.IECDef.numadfn
.. autoattribute:: numadObjects.IECDef.NumGrid
.. autoattribute:: numadObjects.IECDef.numSeeds
.. autoattribute:: numadObjects.IECDef.operatingPoints
.. autoattribute:: numadObjects.IECDef.parDir
.. autoattribute:: numadObjects.IECDef.ratedSpeed
.. autoattribute:: numadObjects.IECDef.sf_fat
.. autoattribute:: numadObjects.IECDef.sf_uts
.. autoattribute:: numadObjects.IECDef.sf_tow
.. autoattribute:: numadObjects.IECDef.SimTime
.. autoattribute:: numadObjects.IECDef.simulinkModel
.. autoattribute:: numadObjects.IECDef.simulinkModelFolder
.. autoattribute:: numadObjects.IECDef.TurbClass
.. autoattribute:: numadObjects.IECDef.ws
.. autoattribute:: numadObjects.IECDef.wd
.. autoattribute:: numadObjects.IECDef.yaw



IEC Methods	
~~~~~~~~~~~~~~~~~~
.. autoattribute:: numadObjects.IECDef.checkInputs
.. autoattribute:: numadObjects.IECDef.setAvgWindSpeed
.. autoattribute:: numadObjects.IECDef.setBladeGageCoordinateRotation
.. autoattribute:: numadObjects.IECDef.setGageLabels
.. autoattribute:: numadObjects.IECDef.setSimFlag
.. autoattribute:: numadObjects.IECDef.runFullLoads
.. autoattribute:: numadObjects.IECDef.setRandomSeeds
   


.. _polarClass:

Polar Class
------------------
.. TODO: properties and methods should be autopopulated, manual for now (below)
.. Kelley: several PolarDef functions are not included within the class, they are appended at the end

.. autoclass:: numadObjects.PolarDef
	:members: 
	:exclude-members: 
	:no-undoc-members: 
	

Polar Properties	
~~~~~~~~~~~~~~~~~~
.. autoattribute:: numadObjects.PolarDef.file
.. autoattribute:: numadObjects.PolarDef.source
.. autoattribute:: numadObjects.PolarDef.titleLine
.. autoattribute:: numadObjects.PolarDef.notes
.. autoattribute:: numadObjects.PolarDef.param
.. autoattribute:: numadObjects.PolarDef.rawlist
.. autoattribute:: numadObjects.PolarDef.rawdata
.. autoattribute:: numadObjects.PolarDef.modopts
.. autoattribute:: numadObjects.PolarDef.modlist
.. autoattribute:: numadObjects.PolarDef.moddata

Polar Methods	
~~~~~~~~~~~~~~~~~~
.. autoattribute:: numadObjects.PolarDef.getRawData
.. autoattribute:: numadObjects.PolarDef.getModData
.. autoattribute:: numadObjects.PolarDef.getParam
.. autoattribute:: numadObjects.PolarDef.matchModList
.. autoattribute:: numadObjects.PolarDef.plotRaw
.. autoattribute:: numadObjects.PolarDef.plotMod
.. autoattribute:: numadObjects.PolarDef.resetModData
.. autoattribute:: numadObjects.PolarDef.updateModData
.. autoattribute:: numadObjects.PolarDef.getModOpts
.. autoattribute:: numadObjects.PolarDef.addModOpts
.. autoattribute:: numadObjects.PolarDef.clearModOpts
.. autoattribute:: numadObjects.PolarDef.apply3DStall
.. autoattribute:: numadObjects.PolarDef.applyExtrap
.. autoattribute:: numadObjects.PolarDef.applyResample
.. autoattribute:: numadObjects.PolarDef.calcDynStall
.. autoattribute:: numadObjects.PolarDef.updateIntrinsic
.. autoattribute:: numadObjects.PolarDef.setIntrinsic
.. autoattribute:: numadObjects.PolarDef.blend
.. autoattribute:: numadObjects.PolarDef.plotInterp
.. autoattribute:: numadObjects.PolarDef.plotIntrinsic
.. autoattribute:: numadObjects.PolarDef.writePolar
.. autoattribute:: numadObjects.PolarDef.cylinderPolar
.. autoattribute:: numadObjects.PolarDef.createModOpts
.. autoattribute:: numadObjects.PolarDef.findZeroCL
.. autoattribute:: numadObjects.PolarDef.findMaxCL
.. autoattribute:: numadObjects.PolarDef.findMaxLoD
.. autoattribute:: numadObjects.PolarDef.findAlphaTrend
.. autoattribute:: numadObjects.PolarDef.recolorplot






