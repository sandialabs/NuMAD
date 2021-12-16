
One-way Coupling to GUI
=======================

.. _NuMAD2p0:


Once the ``BladeDef`` object is defined, it is possible to interface
with the NuMAD GUI. The function that
writes the input file for the GUI from a blade object is:
``source\objects\BladeDef_to_NuMADfile.m``

For example the function can be used as:

.. code-block:: matlabsession
	
    >> BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt')

where ``blade`` is a blade object, ``'numad.nmd'``\ is the desired name
to be given to the NuMAD file and ``'MatDBsi.txt'`` is the desired
material data base name.


