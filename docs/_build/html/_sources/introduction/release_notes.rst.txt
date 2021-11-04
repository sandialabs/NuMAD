.. _intro-release-notes:

Release Notes
=============

.. _NuMADv3:

NuMAD v3.0
----------------
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

-  3D FEA shell analyses have been detached from the GUI and can now be 
   parameterized for various parameter studies and optimization. Data I/O for 
   ANSYS has been automated for mesh generation as well as various analyses, 
   such as tip-deflection, buckling, material rupture, total mass, and
   frequencies. 

-  A technique was developed to determine the design loads. The
   thousands of section forces and moments that occurred during the
   dynamic structural analyses required in the system-level optimization
   were reduced to nine critical load cases.

-  The other newly developed analysis procedure allowed to evaluate
   fatigue damage for every material layer at various cross-sections of
   a blade. Further details on both new analyses are provided in the
   journal article “Part II: 3D FEM design optimization of the rotor
   blades”.

.. Kelley: Link to release on GitHub and add DOI


.. _FutureDev:

Future Development
---------------------

-  Allow for the option to make the blade entirely of solid elements

-  Incorporate adhesive modeling

-  Allow for FEA without the need of commercial FE licenses

-  Probabilistic flaw distributions

-  Incorporate Progressive damage

-  In addition to the ANSYS interface, add capability for a user to use
   other commercial FEA codes such as Abaqus and/or Nastran


Prior Releases
----------------

.. _NuMADv2:

NuMAD v2.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. TODO: Link to release on GitHub and add DOI

* Initial release of NuMAD is on the `NuMAD2p0 static branch <https://github.com/sandialabs/NuMAD/tree/NuMAD2p0>`_ 

* Refer to the former user’s manual in PDF form (`SAND2012-7028 <https://energy.sandia.gov/wp-content/gallery/uploads/NuMAD_UserGuide_SAND2012-7028.pdf>`__).


