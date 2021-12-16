.. _bladeDefAndTerms:

Blade Definition
=================

In NuMAD, a blade is uniquely defined with the ``BladeDef`` object, or blade
object for short. As defined in ``source\objects\BladeDef.m``, many of the properties are parameterized by spanwise location. Refer to
:ref:`bladeClass` for a complete listing of ``BladeDef`` properties.

First and foremost there are *stations*. A station is an airfoil at a
specified span location. The airfoil is partitioned by *keypoints*,
shown in :numref:`bladeKeyPoints`. Various blade properties such as ``blade.leband``,
``blade.teband``, ``blade.sparcapwidth``, and ``blade.sparcapoffset`` help to
position the keypoints precisely. For example, ``blade.leband`` is the
arclength from the *le* keypoint to the keypoint *a*. *Regions* are
defined between the keypoints as listed in :numref:`defineRegions`. An adjacent
station helps define these regions as areas. Spanwise lines emanating
from each keypoint are connected to the corresponding keypoints on an
adjacent station; thus bounding the region with four curves. A suffix of
either HP or LP is added to each region name to distinguish regions on
the high pressure surface verses the low pressure surface. Other airfoil
properties and external blade shape data are defined with the ``AirfoilDef``
class and the ``StationDef`` object respectively.
Usually, the number of stations defined needs to be supplemented for
with interpolated stations.

Material properties, layup information, and thicknesses and widths are
additionally defined in the ``MaterialDef``, ``StackDef``, and ``ComponentDef`` respectively.
Refer to the :ref:`classDefs` for more information.

.. _bladeKeyPoints:
.. figure:: /_static/images/bladeKeyPoints.png
   :width: 6.5in
   :height: 2.23056in

   Relative locations of the blade keypoints.
   
   
.. _defineRegions:
.. table:: Region definition by keypoints (TE-Trailing edge, LE-leading edge)

    +----------------------------------+-----------------------------------+
    | Region Name                      | Bounding Keypoints                |
    +==================================+===================================+
    | LE                               | le & a                            |
    +----------------------------------+-----------------------------------+
    | LE Panel                         | a & b                             |
    +----------------------------------+-----------------------------------+
    | Spar                             | b & c                             |
    +----------------------------------+-----------------------------------+
    | TE Panel                         | c & d                             |
    +----------------------------------+-----------------------------------+
    | TE REINF                         | d & e                             |
    +----------------------------------+-----------------------------------+
    | TE Flatback                      | e & te                            |
    +----------------------------------+-----------------------------------+
    
   
Terminology
--------------

============================ ===================================================
Term or Variable       	 	Definition
============================ ===================================================
HP				High Pressure
LE				Leading Edge
LP				Low Pressure
TE				Trailing Edge
============================ ===================================================







   
