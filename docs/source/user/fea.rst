.. _FEAOps:

Finite Element Analysis Operations
==================================

.. _meshGen:

Mesh Generation
---------------
The function called to generates the FE shell model in ANSYS of a NuMAD
blade is 
 
    >> blade.generateFEA

``source\optimization\structOptimization\layupDesign_ANSYSmesh.m``

.. Note:: 
    It is currently necessary to have created the NuMAD input file from 
    NuMAD 2.0, before attempting to run this function. See :ref:`NuMADv2` for further details.

As an example, the following call would build a mesh for the blade blade
object

.. code-block:: matlabsession

    >> layupDesign_ANSYSmesh(blade,config)

The function also issues commands that calls ANSYS to write a textfile
called ``NLIST.lis``. For each node on the wetted area of the blade, this
file stores the node number and its Cartesian coordinates. This
information is utilized when loads are applied to the FE model in a
later step.

.. _coordinateSystems:

Coordinate Systems
------------------

Loading for the FEA model will typically be derived from FAST or
OpenFAST cross-sectional stress resultant forces and moments at various
spanwise positions. These forces and moments are in coordinate systems
that are aligned with the local principal axes (structural) of the
cross-section in the deformed configuration. Thus, those coordinate
systems change with respect to span. :numref:`numadFASTcsys` contrasts the FAST
results coordinate basis vectors,
:math:`v_{i}^{(j)}\ (j = 1,2,3,\ldots,k)`, with those of the NuMAD loads
system, :math:`w_{i}`, and the ANSYS coordinate system, :math:`x_{i}`,
which are invariant. :math:`v_{i}^{(j)}` is the also known as the blade
coordinate system in FAST/OpenFAST. Here and through the rest of this
document, index notation is used with Latin indices assuming 1, 2, and 3
and Greek indices assuming 1, and 2 (except where explicitly indicated).
Repeated indices are summed over their range (except where explicitly
indicated). 

.. Note::
     The :math:`x_{i}` coordinate system is like the :math:`w_{i}` system 
     but the first axis, :math:`w_{1}` points toward the leading edge 
     instead of the flap direction.

.. _numadFASTcsys:
.. figure:: /_static/images/numadFASTcsys.png
   :width: 5.85771in
   :height: 4.10039in
   
   Comparison of the FAST results coordinate system with that of the NuMAD and ANSYS coordinate systems.

A small angle is assumed to exist between
:math:`v_{3}^{(j)}\ (j = 1,2,3,\ldots,k)` and :math:`w_{3}` and is
unaccounted for in NuMAD. Thus, it is assumed that
:math:`v_{3}^{(j)} = w_{3} = x_{3}`. Thus the :math:`v_{i}^{(j)}` and
the :math:`w_{i}` are related by

.. math:: w_{i} = {C_{iq}v}_{q}^{(j)}

where :math:`C_{iq}` are the direction cosines which can be
stored in a matrix as shown below
        
.. math:: C = \begin{bmatrix}
              \text{cos}(\mu) & \text{sin}(\mu) & 0 \\
              -\text{sin}(\mu) & \text{cos}(\mu) & 0 \\
              0 & 0 & 1
              \end{bmatrix}   
    :label: DCM

      
where :math:`\mu` is the so-called structural twist and is illustrated
in :numref:`FASTvsNuMADcsys`.

.. _FASTvsNuMADcsys:
.. figure:: /_static/images/FASTvsNuMADcsys.png
   :width: 3.15956in
   :height: 1.77725in
   
   Time history moment data from the :math:`v_{\alpha}` coordinate system to the :math:`w_{\alpha}` intermediate coordinate system.

These angles are obtained by `PreComp <https://www.nrel.gov/wind/nwtc/precomp.html>`__. This coordinate transformation is
implementation in ``loadFASTOutDataGageRot``.

Before transforming the loads data to the :math:`x_{i}` system, various
data processing operations occur in intermediate coordinate systems as
described in the next section.

.. _analDir:

Analysis Directions
-------------------

Increasing the fidelity from a beam model in FAST/OpenFAST to a shell
model also increases the computational cost. Thus, it would be
cumbersome to conduct a time history structural analyses with a shell
model for each DLC; it would be impracticable on each iteration of an
optimization loop. Therefore, it is necessary to construct a reduced set
of load cases that is representative of the most critical loads in the
time history analysis.

Two load cases were deemed necessary to properly analyze the most
critical loads; those needed to evaluate the maximum tip deflection and
those needed for evaluating blade failure. Here, blade failure consists
of ultimate failure, buckling, and fatigue failures. Both cases utilize
the FAST/OpenFAST resultant axial forces :math:`F_{3}^{w}`, both bending
moment components :math:`M_{\alpha}^{w}`, and the twisting moment,
:math:`M_{3}^{w}`. The superscript indicates the reference frame. These
will be referred to collectively with :math:`P_{i}\ (i = 1,2,\ldots,4)`,
a generalized load vector, where :math:`P_{i} = \begin{bmatrix}
F_{3}^{w} & M_{1}^{w} & M_{2}^{w} & M_{3}^{w} \\
\end{bmatrix}^{T}`. 

.. Note:: 
    Note that :math:`P_{i}` varies with span-wise location and that resultant shear forces were assumed to be negligible for establishing load equivalency for both load cases. 

For the tip deflection case, the time at which the maximum tip deflection from the
beam model, :math:`t^{*}\ `, was determined. Then the :math:`P_{i}`, at
a given cross section was defined by

.. math::

   P_{i} = \begin{bmatrix}
   F_{3}^{w}(t^{*}) & M_{1}^{w}(t^{*}) & M_{2}^{w}(t^{*}) & M_{3}^{w}(t^{*}) \\
   \end{bmatrix}^{T}

The components were then transformed to the :math:`x_{i}` system.

Since IEC allows for lower factors of safety if numerous analysis
directions are considered, :math:`n_{\theta}` analysis directions can be
considered by letting params.momentMaxRotation\ :math:`{= n}_{\theta}`
in the ``runIEC_ipt.m`` file.

Thus the loads used to evaluate blade failure are obtained by letting
:math:`\theta`, as defined in :numref:`loadDirections`, vary from 0 to 360 deg. in
increments of :math:`n_{\theta}/360`; yielding :math:`n_{\theta}` FE
load cases. The results for the load directions
0\ :math:`\leq \theta`\ <180 deg. were obtained by

:math:`P_{i} = \begin{bmatrix}
\text{max}(F_{3}^{w}(t)) \\
\text{max}(M_{1}^{y}(t))\cos(\theta) \\
\text{max}(M_{1}^{y}(t))\sin(\theta) \\
\text{max}(M_{3}^{w}(t)) \\
\end{bmatrix}` 0\ :math:`\leq \theta`\ <180

where

.. math:: M_{1}^{y} = M_{1}^{w}\left( t,w_{3} \right)\cos(\theta) + M_{2}^{w}(t,w_{3})\sin(\theta)

Load directions from 180\ :math:`\leq \theta`\ <360 deg. are, however,
obtained by

:math:`P_{i} = \begin{bmatrix}
\text{max}(F_{3}^{w}(t)) \\
\text{min}(M_{1}^{y}(t))\cos(\theta) \\
\text{min}(M_{1}^{y}(t))\sin(\theta) \\
\text{max}(M_{3}^{w}(t)) \\
\end{bmatrix}` 0\ :math:`\leq \theta`\ <180

.. Note:: 
    In the above equation, the minimum of :math:`M_{1}^{y}` is found 
    instead of the maximum but :math:`\theta` still ranges from 0 to 180 deg. 
    Unlike the resultants used for the deflection analysis which all occurred 
    at a specific time in the OpenFAST simulations, each of the axial force, 
    torsion, and bending moment resultants along the span could possibly come
    from different times. Thus, it is an artificial distribution suitable for
    design use. Also unlike the deflection case, the two bending-moment 
    components (i.e. flap-wise and edge-wise bending) at a span location
    were projected onto 12 directions as defined by :math:`y_{1}` in 
    :numref:`loadDirections`. All of the :math:`P_{i}\ `\ then transformed 
    to the :math:`x_{i}` system in ``beamForceToAnsysShell``.

.. _loadDirections:
.. figure:: /_static/images/loadDirections.png
   :width: 3.06468in
   :height: 1.72388in
   
   The definition of the :math:`y_1` axis and the angle :math:`\mathbf{\theta}`.

.. _applyLoads:

Loads Application
-----------------

Blade loads can either come from the cross-sectional stress resultant
time-histories from FAST or user defined. Two load types are supported
for loads coming from a FAST; loads at a particular time and extremum
loads used to evaluate the limit states of the blade. Both account for
all three cross-sectional resultant moments along the span as well as
the resultant axial forces. However, the resultant transverse shear
forces acting on the blade in the FAST model are not transferred the
ANSYS model. Moreover, both create the ``loadsTable`` variable needed by
``layupDesign_ANSYSanalysis``, the main FEA script described is subsequent
section. ``getForceDistributionAtTime.m`` handles can be called for the
loads at a given time while ``FastLoads4ansys.m``} builds the ``loadsTable`` for
each analysis direction. The ``loadsTable`` variable
is\ :math:`\ 1 \times n` MATLAB cell array where :math:`n` is
the number of analysis directions. 

.. Note:: 
   ``getForceDistributionAtTime.m`` will give a :math:`1 \times 1` cell array since the loading is the *actual* loading at a particular time.). 

A user may apply loads other than the loads at particular time or the extremum loads with a user
defined ``loadsTable``. This can be done by creating the structure of the
``loadsTable`` variable manually.

Example use for ``FastLoads4ansys.m`` is as follows:

.. literalinclude:: loadsExample.m
   :language: matlab

For bending moments on a blade, ``getForceDistributionAtTime.m`` and
``getForceDistributionAtTime.m`` convert the moment distributions to
transverse force distributions. :numref:`momentsToForces` shows the known spanwise
bending moments, :math:`M_{i}` acting at a distance :math:`z_{i}` from
the blade root and the forces to be solved, :math:`F_{i}` acting at
:math:`{\overline{z}}_{i}`, where
:math:`{\overline{z}}_{i} = 1/2(z_{i} + z_{i + 1})`.

.. _momentsToForces:
.. figure:: /_static/images/momentsToForces.png
   :width: 6.5in
   :height: 3.65625in

   Freebody diagram used to determine transverse loads from a given distribution of moments along the blade span.


From static equilibrium, the forces to be applied are found by solving
the following linear system of equations for
:math:`F_{i}\ (i = 1,2,3,\ldots,k)`

.. math:: M_{i} = \sum_{j = i}^{k}{F_{j}({\overline{z}}_{j} - z_{i})}

where :math:`k` is the number of cross-sections with resultant moment
data.

The load distributions are then transferred to nodal loads by
``beamForceToAnsysShell.m``. ``beamForceToAnsysShell.m`` generates a text file (usually called forces.src)
containing the APDL commands to apply nodal forces for every node on the
wetted area of the blade. Details of the approach are found in Ref. [1]
but with modifications to add axial loads. 

.. Note:: 
    The script assumes that the input forces are in the FAST coordinate system, 
    :math:`w_{i}`, so a conversion is performed in the script to the 
    ANSYS coordinate system, :math:`x_{i}`.

A vector, :math:`G_{i}`, such as a force or a
moment, in the :math:`w_{i}` coordinate system is transformed to the
:math:`x_{i}` with

.. math:: g_{i} = C_{ij}G_{j}

where :math:`C_{ij}` is defined in Eq. :eq:`DCM` but with setting
:math:`\mu = - 90\ `\ deg. Example usage for the first ``loadsTable`` load
case is as follows

.. literalinclude:: loadsTableExample.m
   :language: matlab

where maptype is either ``'map2D_fxM0'`` or ``'map3D_fxM0'``. Note that these
steps are not usually required with using ``layupDesign_ANSYSanalysis``, the
main analysis script. These instructions are for building the ``loadsTable`` variable
and applying nodal forces to an FE model in a stand-alone manner.

.. _analysisScript:

Analysis Script
---------------

``layupDesign_ANSYSanalysis.m`` is the main analysis script. A user may call
on it to

-  obtain load factors from global instabilities due to linear or
   nonlinear buckling

-  obtain load factors from local instabilities due to face-sheet
   wrinkling

-  find the maximum value of a failure index

-  stress data required for fatigue analyses

-  cross sectional force and moment results vs spanwise position

-  average cross-sectional rotations and displacements along the span

-  strain along the spar caps vs. spanwise position

-  blade mass

The required inputs are a blade object, the ``loadsTable``, and the ``config``
variable. The script knows which analysis types it should run based on
config. See :numref:`configTable` for the structure and usage of the ``config``
variable.


.. _configTable:
.. csv-table:: Structure and usage of the ``config`` variable.
   :file: configTable.csv
   :widths: 1, 1, 3
   :header-rows: 1
   

``meshFile`` is the name of the ANSYS model database to begin the analysis
from. It defaults to ``'master'``. analysisFileName is the name of the new
ANSYS model database which will store the analysis results. It defaults
to ``'FEmodel'``. np is the number of processors (CPUs) ANSYS should use for
each solve. It defaults to 1.

The output and it is a variable length struct. Depending on
which analysis flags are active, results can be accessed with

.. code-block:: matlabsession

    >> result=layupDesign_ANSYSmesh(blade,config)
    >> result.globalBuckling
    >> result.localBuckling
    >> result.deflection
    >> result.failure
    >> result.fatigue
    >> result.resultantVSspan
    >> result.mass


.. _linearFEA:

Linear-Static Analysis
~~~~~~~~~~~~~~~~~~~~~~

A linear-static analysis is defined and solved by default. Since
linear-buckling is a notable advantage of using a shell model for the
blade, prestress is activated in this load step in preparation for a
linear buckling analysis. The effects of inertia are counterbalanced
with inertia relief.

Depending on the ``config`` flags, a number analysis outputs can be
examined.

.. _deflectionFEA:

Deflection
~~~~~~~~~~

Setting ``config.ansys.deflection = 1`` will cause ``result.deflection`` to
hold the average deflection of each cross-section at every ``blade.ispan``
location. For example, access results from the second ``loadsTable`` case as

``result.deflection{2}``

The result is a matrix that has as many rows are ``blade.ispan`` and 6
columns. The first three columns are average cross-sectional
displacements in the :math:`x_{i}` system. Columns 1, 2, and, 3,
correspond to :math:`{\overline{u}}_{1}`, :math:`{\overline{u}}_{2}`,
and :math:`{\overline{u}}_{3}` respectively. Columns 4, 5, and 6
correspond to :math:`{\overline{\theta}}_{1}`,
:math:`{\overline{\theta}}_{2}`, and :math:`{\overline{\theta}}_{3}`,
respectively.

.. _failureFEA:

Failure
~~~~~~~

If material failure is being considered in an analysis, set
``config.ansys.failure`` equal to any one of the entries in :numref:`failureCriteriaOptions`. An
empty string will skip any material failure analyses.


.. _failureCriteriaOptions:
.. csv-table:: Supported failure criteria and associated ``config`` flags.
   :file: failureCriteriaOptions.csv
   :widths: 2, 2
   :header-rows: 1



Depending on the criterion used, it is necessary for the ``materialDef`` to
have the appropriate material properties defined in ``blade.materials``.

Since each element section is a layered composite, by default,
``result.failure`` will hold the maximum failure index of all elements and
every layer.

.. _globalBucklingFEA:

Global Buckling
~~~~~~~~~~~~~~~

A nonzero number for ``config.ansys.analysisFlags.globalBuckling`` will call
ANSYS to run a buckling analysis. The buckling analysis can either be
linear or nonlinear. The nonlinear case will be activated by creating
``config.ansys.analysisFlags.imperfection`` and setting it to a nonempty
number array. ANSYS will then be directed to perform an eigenvalue
buckling analysis. The load which results in nonconvergence will be
reported as the nonlinear buckling load factor.

For the nonlinear case, it is necessary to introduce a synthetic
imperfection for the analysis to provide reasonable results. The
imperfection introduced is of the geometric kind and corresponds to the
buckled mode shapes obtained in the eigenvalue buckling analysis.
Therefore, an eigenvalue buckling analysis precedes the nonlinear static
analyses. Currently, for nonlinear buckling, it is assumed here that
buckling will be in the aeroshell (not the web). Thus, each mode shape
is scaled such that the maximum displacement in the :math:`x_{2}`
direction (flapwise) is equal to config.ansys.analysisFlags.imperfection
in value. Robustness could be increased if the script can select whether
to find the maximum displacement in the :math:`x_{2}` (buckling in
aeroshell) or :math:`x_{1}` direction (buckling in the webs).

Finally, it is supposed that the mode shape corresponding to the lowest
load factor from the eigenvalue analysis will not always cause the
lowest load factor from the nonlinear case. Thus, a nonlinear static
analysis for each requested mode is performed and the minimum load
factor is reported.

.. _localBucklingFEA:

Local Buckling 
~~~~~~~~~~~~~~

Sandwich panels typically consist of a relatively soft core layer
*sandwiched* between two faces sheets. For the parts of blade that are
described as sandwhich panels, local instabilities can be checked in
NuMAD with a strip theory model from Ref. [2]. In particular, this check
examines if the outermost facesheet wrinkles under compression.

To activate this check, set ``config.ansys.analysisFlags.localBuckling``
equal to the name of the core material in your BladeDef as a string
(e.g. ``config.ansys.analysisFlags.localBuckling = 'Balsa'``). Otherwise an
empty string will skip the analysis. Currently, all ply angles are
limited to zero (i.e. no off-axis layups).


.. _fatigueFEA:

Fatigue
~~~~~~~

If ``config.ansys.fatigue`` is a nonempty string, then ANSYS will be
directed to write model data and field output to text files for fatigue
post-processing in MATLAB. Namely, ``Elements.txt``, ``Materials.txt``,
``Strengths.txt``, ``Sections.txt``, and plate strain measures. The results from running
the fatigue post-processor will be a :math:`1 \times n_{\theta}/2` cell
array, where :math:`n_{\theta}` is the number of cells in the ``loadsTable``
variable (where each cell holds the loads for a single analysis
direction). This is because stress results from two orthogonal loading
directions are utilized for a single evaluation of fatigue damage (e.g.
one fatigue evaluation is the combined effect of flap loads and edge
loads). Furthermore, it assumes that the ``loadsTable`` is arranged in
ascending order for the loads direction angle. This is already accounted
for if ``loadsTable`` was constructed from ``FastLoads4ansys.m``.
Make sure to set ``params.fatigueCriterion = 'Shifted Goodman'`` 
in the ``runIEC_ipt.m`` file.

Fatigue damage is computed at the root and at the locations of the FAST
gages. In this document, these will be referred to as the *spanwise
fatigue locations*. Fatigue results are restricted to the spanwise
fatigue locations because the fatigue analysis relies on Markov matrices
and these matrices need cycle counts from time-history data. The time
histories that are available from FAST are at the root and the FAST
gages.

The end-user has the ability to either obtain a single fatigue damage
fraction at each spanwise fatigue location or several fatigue damage
fractions each corresponding to a region in :numref:`bladeKeyPoints`. The former is
accomplished by setting ``config.ansys.analysisFlags.fatigue = ["ALL"]``.
The latter is accomplished by setting ``config.ansys.analysisFlags.fatigue``
to a MATLAB string array where any and all Region Names in :numref:`regionNames` are allowed.


.. _regionNames:
.. csv-table:: Region Names to specify in config for fatigue analyses.
   :file: regionNames.csv
   :widths: 2,5
   :header-rows: 1


For example, suppose that the fatigue damage was only desired at the
spar caps and the reinforcement locations one would set

.. code-block:: matlabsession

    >> config.ansys.analysisFlags.fatigue = ["HP_TE_FLAT","HP_TE_REINF","HP_SPAR","HP_LE","LP_LE","LP_SPAR","LP_TE_REINF","LP_TE_FLAT"]

If ``"ALL"`` is included along with other regions, the other regions are
ignored and the analysis carries forth as it would for ``"ALL"``.

The result will be a :math:`1 \times n_{\theta}/2` cell array. Each cell
will hold a struct with the following field names: ``fatigueDamage``,
``criticalElement``, ``criticalLayerNo``, and ``criticalMatNo``. Each of these field
names allow access to a :math:`n_{\text{region}} \times 10` matrix where
:math:`n_{\text{region}}` is the length of
``config.ansys.analysisFlags.fatigue``. The results are organized such that
the :math:`i_{\text{th}}` row number corresponds to fatigue damage
results of the :math:`i_{\text{th}}` Region Name from
``config.ansys.analysisFlags.fatigue``. The columns correspond to the
fatigue damage spanwise locations where column 1 is the results at the
root. Subsequent columns correspond to successively farther gage
locations.

The results from the first fatigue evaluation can visualized as shown in
:numref:`fatigueDamageFractionExample`.

.. _fatigueDamageFractionExample:
.. figure:: /_static/images/fatigueDamageFractionExample.bmp
   :alt: Diagram Description automatically generated
   :width: 4.17467in
   :height: 3.131in

   Fatigue damage fraction visualization.
	

.. _localFieldsFEA:

Local Fields
~~~~~~~~~~~~

If ``config.ansys.analysisFlags.localFeilds`` is activated, then the plate
strains and curvatures for every element will be written to
plateStrains-all-1.txt.

After layupDesign_ANSYSanalysis completes execution, access the local
fields with

.. literalinclude:: localFieldsExample.m
   :language: matlab

``elNo`` is an scalar or array of element numbers for
which to extract the local fields. Set ``coordSys='local'`` to obtain
results in the local layer coordinate system or set ``coordSys='global'`` to
obtain results in the element coordinate system.

For example if one sets ``elNo=[1950,2558]`` then, ``myresult`` is a struct with
two fields:

.. code-block:: matlabsession

    >> myresult
	     element1950: [1x1 struct]
	     element2558: [1x1 struct]

Looking at the first struct shows:

.. code-block:: matlabsession

    >> myresult.element1950
         x3: [8x1 double] 
         eps11: [8x1 double] 		 
         eps22: [8x1 double]		 
         eps33: [8x1 double]		 
         eps23: [8x1 double]		 
         eps13: [8x1 double]		 
         eps12: [8x1 double]		 
         sig11: [8x1 double]		 
         sig22: [8x1 double]		 
         sig12: [8x1 double]		 
         matNumber: [8x1 double]

where ``x3`` is the thickness coordinate, ``eps`` are the strains and sig are the
stresses. 

.. Note:: 
    Note the ``x3`` value that is zero corresponds to the reference surface in the shell model. 
	
In this example the shell offset in ANSYS was set to ``BOT`` so the zero appears the end:

.. _localFieldsDataExplanationFigure:
.. figure:: /_static/images/localFieldsDataExplanationFigure.png
   :width: 5.03004in
   :height: 2.30597in

   Explanation of how data for local fields is ordered.


A plotting utility has been made to visualize the local fields. This can
be achieved for example:

.. literalinclude:: localFieldsPlotExample.m
   :language: matlab

.. _localFieldsPlotExample:
.. figure:: /_static/images/localFieldsPlotExample.png
   :width: 6.5in
   :height: 2.98472in

   Local fields example plots.
	

.. _resultantsFEA:

Resultant Forces and Moments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The resultant forces and moments due to stresses on cross-sections along
the span can be obtained by activating
``config.ansys.analysisFlags.resultantVSspan``. ``result. resultantVSspan{i}``
will yield a :math:`n \times 7` matrix, where :math:`n` is approximately
equal to the blade length divided by blade.mesh. The first column is the
distance along :math:`x_{3}` that locates the cross-section. The next
three columns correspond to resultant forces in the :math:`x_{1}`,
:math:`x_{2}`, and :math:`x_{3}` directions (in that order). The next
three columns correspond to resultant moments about in the
:math:`x_{1}`, :math:`x_{2}`, and :math:`x_{3}` directions (in that
order). For example, access results such as

.. code-block:: matlabsession

    >> result.resultantVSspan{i}

where i is the i\ :sup:`th` load case in ``loadsTable``.
