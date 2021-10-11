.. _intro-release-notes:

Release Notes
=============

Current Release
----------------

NuMAD v3.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The current release now incorporates structural optimization, associated
structural analyses, and the move to object-oriented data structures.
The exclusive use of the GUI in prior versions prevented automation in
optimizations. Thus, moving to these data structures enabled
optimization. Since the GUI functionality is advantageous in certain
situations, its functionality is still retained. Current documentation primarily
describes usage of the structural optimization, associated structural
analyses, and the object-oriented data structures. 

-  A new capability also includes the ability to accept input from the
   International Energy Agency (IEA) Wind Task 37 blade ontology. This
   is known to have enhanced collaboration.

-  Although analysis of tip-deflection, buckling, material rupture,
   total mass, and frequencies were possible in NuMAD 2.0, the new
   version automates the transfer of those data from ANSYS back to NuMAD
   for the optimizer or post-processing.

-  A technique was developed to determine the design loads. The
   thousands of section forces and moments that occurred during the
   dynamic structural analyses required in the system-level optimization
   were reduced to nine critical load cases.

-  The other newly developed analysis procedure allowed to evaluate
   fatigue damage for every material layer at various cross-sections of
   a blade. Further details on both new analyses are provided in the
   journal article “Part II: 3D FEM design optimization of the rotor
   blades”.

NuMAD v2.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Initial release of NuMAD is on the NuMAD2p0 static branch (originally on `Sandia National Laboratories (Sandia) <https://energy.sandia.gov/programs/renewable-energy/wind-power/>`, now available on GitHub)

* Refer to the former user’s manual in PDF form (`SAND2012-7028 <https://energy.sandia.gov/wp-content/gallery/uploads/NuMAD_UserGuide_SAND2012-7028.pdf>`__).



