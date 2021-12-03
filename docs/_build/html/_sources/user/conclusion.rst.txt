.. _troubleshooting:

Troubleshooting
===============

NuMAD will repeatedly print a message in the MATLAB command window to the effect of:

.. code-block:: matlabsession

	``Waiting for ANSYS to <`do something`>``
	``Waiting for ANSYS to <`do something`>``
	``Waiting for ANSYS to <`do something`>``

Common things to check are the ANSYS `.err` and `.log` for clues. 

.. _KnownIssues:

Known Issues
============

Trailing Edge Issues

-  There is a need to allow for the existence of a nonzero trailing edge
   thickness. The only way to achieve this in the current release is to
   specify that the airfoil in question is a flatback.
   
-  Trailing edge overlap

.. _TEoverlap:
.. figure:: /_static/images/TEoverlap.png
   :width: 5.85771in
   :height: 4.10039in

-  Often the element quality is poor at various parts of the trailing
   edge. Possible solutions could be:

   -  Alter airfoil geometry when reading in to make sharp airfoils for
      small trailing edge thicknesses

   -  Incorporate solid elements to represent trailing edge adhesive

Other Issues

-  The placement of the spar cap must never exceed the bounds of the
   spar cap keypoints.

-  ``blade.updateBOM`` fails for an unknown reason if the layer quantity
   increases rather than decreases toward the tip.



.. _References:

References
==========

1. Berg, Jonathan, Joshua Paquette, and Brian Resor. "Mapping of 1D beam
   loads to the 3D wind blade for buckling analysis." *52nd
   AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics and Materials
   Conference 19th AIAA/ASME/AHS Adaptive Structures Conference 13t*.
   2011.

2. Fagerberg, Linus, and Dan Zenkert. "Effects of anisotropy and
   multiaxial loading on the wrinkling of sandwich panels." *Journal of
   Sandwich Structures & Materials* 7.3 (2005): 177-194.

3. Berg, Jonathan C., and Brian R. Resor. "Numerical manufacturing and
   design tool (NuMAD V2. 0) for wind turbine blades: user's guide."
   *Sandia National Laboratories Technical Report, SAND2012-7028*
   (2012). 


.. TODO: use userGuide.bib publications to reference citations
   