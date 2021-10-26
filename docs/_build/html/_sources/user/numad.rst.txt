
Use Cases
=====================

.. _NuMAD2p0:

Coupling to NuMAD 2.0
----------------------
Once the ``BladeDef`` object is defined, it is possible to interface
with prior versions of NuMAD which were GUI-centric. The function that
writes the input file to older versions of NuMAD from a blade object is:
``source\objects\BladeDef_to_NuMADfile.m``

For example the function can be used as:

.. code-block:: matlabsession
	
	``BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt')``

where ``blade`` is a blade object, ``'numad.nmd'``\ is the desired name
to be given to the NuMAD file and ``'MatDBsi.txt'`` is the desired
material data base name.


.. _GUI:

NuMAD GUI Mode
-----------------
In this version of NuMAD, the graphical user interface (GUI) can still
be accessed the same as it was in prior releases. Refer to the former user manual
(`SAND2012-7028 <https://energy.sandia.gov/wp-content/gallery/uploads/NuMAD_UserGuide_SAND2012-7028.pdf>`__) for
detailed instructions on how to use the GUI . For operating NuMAD exclusivley with the GUI refer
to the :ref:`intro-release-notes` on NuMAD v2.0.