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

-  3D FEA shell analyses have been partially detached from the GUI and can now be 
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

-  Click `here <https://github.com/sandialabs/NuMAD/releases/tag/v3.0>`_ for NuMAD v3.0
.. Kelley: add DOI


.. _FutureDev:

Possible Future Development
---------------------------

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
.. TODO: add DOI


NuMAD first became available in the early 2000’s. It was originally written in Tcl, chosen for its ability to easily manipulate graphics and port information in and out of the main code. NuMAD has always been designed to output files for ANSYS. The latest release of the previous NuMAD code is dated 31 March 2010.

Beginning with the last release of the previous NuMAD, its authors realized that blade design was quickly heading in a direction where a computationally intense, and platform independent programming environment was needed for NuMAD. Mathworks MATLAB® was chosen as the environment because of its widespread use, its computational power, its flexibility in operating systems, its graphical capabilities, its popularity with researchers and students, and its ability to compile as well as to run from raw source code. The current release of NuMAD also includes some advanced capabilities which go beyond the basic connection to creation of the finite element model: tabularized input format, swept blades, blades with pre-bend, output for CFD mesh creation and output of blade cross section properties.

NuMAD v2.0 is available in two forms: 1) compiled MCR executable and 2) raw MATLAB source files. The MCR executables are in binaries.zip. Note that the MCR executables were too large to be placed in the "bin" folder of the source code. They as well as the rest of the "bin" files are in the binaries.zip asset. Thus, you can ignore the files in "bin" folder of the source code. Users are encouraged to study the raw source files to understand how NuMAD works. When there are capabilities that presently do not exist, the authors encourage users to work together with the authors to write modules that accomplish the tasks. It is the hope of the authors that the accessibility of the MATLAB source code will, in the long run, enable an extremely capable tool that will continue to enable DOE’s goals and benefit the entire wind industry.

* Initial release of NuMAD in MATLAB: `Click here <https://github.com/sandialabs/NuMAD/releases/tag/v2.0>`_ 

* Refer to the former user’s manual in PDF form (`SAND2012-7028 <https://energy.sandia.gov/wp-content/gallery/uploads/NuMAD_UserGuide_SAND2012-7028.pdf>`__).

