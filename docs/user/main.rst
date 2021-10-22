
.. _bladeDefAndTerms:

Blade definition and terminology
================================

In NuMAD, a blade is uniquely defined with the ``BladeDef`` object or blade
object for short. As defined in

``source\preNuMAD\BladeDef.m``,

many of the properties are parameterized by spanwise location. Refer to
:numref:`bladeDefTable` in the :ref:`appendix` for a complete listing of ``BladeDef`` properties.

First and foremost there are *stations*. A station is an airfoil at a
specified span location. The airfoil is partitioned by *keypoints* as
shown in :numref:`bladeKeyPoints`. Various blade properties such as ``blade.leband``,
``blade.teband``, ``blade.sparcapwidth``, and ``blade.sparcapoffset`` help to
position the keypoints precisely. For example, blade.leband is the
arclength from the “le” keypoint to the keypoint “a”. *Regions* are
defined between the keypoints as listed in :numref:`defineRegions`. An adjacent
station helps define these regions as areas. Spanwise lines emanating
from each keypoint are connected to the corresponding keypoints on an
adjacent station; thus bounding the region with four curves. A suffix of
either HP or LP is added to each region name to distinguish regions on
the high pressure surface verses the low pressure surface. Other airfoil
properties and external blade shape data are defined with the ``AirfoilDef``
class and the ``StationDef`` object respectively.


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

Usually, the number of stations defined needs to be supplemented for
with interpolated stations.

Material properties, layup information, and thicknesses and widths are
additionally defined in the

``MaterialDef``, ``StackDef``, and ``ComponentDef`` respectively.

.. _bladeKeyPoints:
.. figure:: /_static/images/bladeKeyPoints.png
   :width: 6.5in
   :height: 2.23056in

   Relative locations of the blade keypoints.

.. _gettingStarted:

Getting Started with NuMAD
==========================

The code containing all open-source tools for NuMAD can be downloaded or
cloned as a repository from the `NuMAD
GitHub <https://github.com/sandialabs/NuMAD>`__. There
are 3 main folders within the root ``source`` directory, named ``NuMAD_toolbox``,
``preNuMAD``, and ``rotor_optimization``. It is advisable to save these
in a directory called ``DesignCodes`` on, for example, the C drive if possible. If another
path location is desired, some path definitions in the ``runIEC_ipt.m``
input settings script may need to be modified, as described further in the
:ref:`AeroelasticSimRunIEC` section, and the :ref:`appendix`. Before beginning any design or analysis of a
blade, the paths of these folders should be added to the MATLAB domain
of working directories, using the script ``paths.m``. All tools and
functions available in NuMAD 3.0 will then be available to invoke from a
script or the command window. If you intend to use NuMAD in conjunction with ANSYS then the global variable 
``ANSYS_Path`` needs to be modified to your ANSYS working directory file path name.

The ``NuMAD_toolbox`` folder contains basic functions and operations needed
for performing analysis with packages such as PreComp, BPE, and ANSYS.
The ``preNuMAD`` folder mainly contains the class definition of the blade
object, which stores the geometric, airfoil and material data for a
given blade design. The ``rotor_optimization`` folder contains bundles of
tools for several purposes, including ``runIEC``, as explained further in the
:ref:`AeroelasticSimRunIEC` section, file processing functions for input and output from other
programs such as FAST and Crunch, and setup and execution of ANSYS
models for analyzing quantities like material rupture, fatigue, and
buckling under loading. 

The following sections go through the basics of putting together and
analyzing a blade design in NuMAD, from creating a blade object, to
modifying its parameters and attributes, to generating structural models
for analysis. Tips for implementation are given to ensure the smoothest
execution.

.. _bladeGen:

Blade Generation
----------------

There are several ways to go about constructing a blade model in NuMAD
3.0. In many ways the most intuitive approach is by using the graphical
user interface carried over from NuMAD 2.0 (for detailed instructions,
please refer to the NuMAD 2.0 user manual [3]). Although the current
version is designed not to be reliant on this GUI, it is still supported
and can be a useful tool if designing a blade model from scratch.

A more automated way to generate a blade model is by reading a ``.yaml``
file, and using the data stored within to populate the model definition.
The ``.yaml`` format contains all the geometric, aerodynamic and material
information needed to define a blade structure, and is widely used in
the field, making it a convenient choice for the source data file.
Several examples of ``.yaml`` files are given in the ``examples`` directory
of the GitHub repository.

To read a ``.yaml`` file’s data into NuMAD from a MATLAB script or command
window, first save the ``.yaml`` file in the working directory for the blade
model. Typically, each blade design should have its own working
directory, which contains aerodynamic/airfoil data, FAST input files
(``*.fst`` extension) and subfolders for NuMAD operations and FAST
simulation output files (``*.out`` extension). Once the ``.yaml`` file is saved
alongside these items, it can be read by creating a new blade object of
type ``BladeDef`` and calling the reader function, as shown:

    >> ``blade = BladeDef``

    >> ``blade.readYAML(<yamlfilename>)``

After entering these commands, the blade model data will be stored in
the object called ``blade``, which can be printed out and modified as
desired as further explained in :ref:`bladeGen`. Many types of analysis and
operations can be performed once a blade model is stored as a blade
object in this way. If it is desired to run aeroelastic simulation with
``runIEC``, or to generate loads from FAST-generated output files as
described in :ref:`FEAOps`, it is advisable to update the FAST input files
in the blade’s main working directory to ensure that certain quantities
within are consistent with the data stored in the ``.yaml`` file, such as
prebend, presweep, and structural twist. This can be done using the
commands

``runIEC_ipt``

``updateFASTFromBLADEDef(params,blade)``

The ``runIEC_ipt.m`` script initializes a data structure named ``params``, which
stores variables related to the analysis and simulation (``\examples\runIEC_ipt--EXAMPLE.m``), and is passed along with the blade object into the update
function.

After all desired analyses and modifications are complete, a new,
updated ``.yaml`` file can be generated to represent the optimized, or
re-designed blade. To do this issue the command:

``writeYAML(blade,<newyamlfilename>);``

The new file name should be of the form ``<originalfile>_mod.yaml``, adding
an ‘_mod’ extension to the name of the original ``.yaml`` file from which
the blade model was generated. The new file will be written into the
blade model’s working directory alongside the original.

Generating blades from ``.yaml`` files is useful for streamlining analysis
and optimization processes, since all operations can be called from a
MATLAB script, without depending on the graphical user interface or any
manual input during execution.

.. _bladeVarAssign:

Blade Variable Assignment
-------------------------

Any blade design process involves setting and modifying characteristics
of the blade’s geometric, structural, and material properties in some
way. The NuMAD blade object contains a collection of variables that
represent these properties, and can be set and modified by value
assignment within a MATLAB script or in the command line. A
comprehensive list of all public variables in the BladeDef class used in
NuMAD is given in the :ref:`appendix`.

Here we give some highlighted examples of key variables within the blade
object and their basic access. The overall shape of the blade is defined
largely by the stations (access: ``blade.stations``). Each station
contains several variables whose values can be edited. For instance, if
a blade model was generated by reading a ``.yaml`` file as described in
:ref:`bladeGen`, and it was desired to set the spanwise position of the
second station at 3.5 meters, the command

``blade.stations(2).spanlocation = 3.5``

could be used. Many variables are arrays with multiple values, and can
be set according using standard MATLAB syntax. The coordinates of the
points defining the outer airfoil shape at a given station, for example,
are stored in the airfoil object at each individual station as an :math:`N X 2`
array, and can be set as follows:

``blade.stations(2).airfoil.coordinates = [X1, Y1; …``

``X2, Y2; …``

``…``

``XN, YN]``

There are several properties that each define some aspect of the blade’s
shape with a value at any given spanwise location, including chord
length, angle of twist, aerodynamic center, sweep and prebend. These can
be set at any number of spanwise points, with the variable ‘span’
specifying their locations. If a user wanted to, say, set the prebend of
the blade to some constant :math:`k` times cube of the spanwise location,
specified at 10 equally spaced points, they could set

``blade.span = linspace(0,<bladeLength>,10)’;``

``blade.prebend = k*blade.span.^3;``

The bulk of the structural properties of the blade’s components are
stored in blade.components variable. A single component contains a name,
a material ID number, labels representing the points it spans between
according to :numref:`bladeKeyPoints`, and a control point array, called ``cp``. The
control point array specifies the thickness of the given component at
every spanwise location, expressed in number of layers (the actual
thickness of a layer is defined by the material object it corresponds
to, shown shortly). Suppose component 3 in the blade was the suction
side spar cap, and it was desired to vary the thickness linearly from 10
layers at the root to 2 layers at the tip, say 50 meters span. The user
could set

``blade.components(3).cp = [0, 10; …``

``50, 2];``

The width of the spar caps and the leading edge and trailing edge bands
are single nominal values for the entire length of the blade, stored in
the variables ``blade.sparcapwidth``, ``blade.leband`` and ``blade.teband``
respectively.

The data defining the properties of all the materials used throughout
the blade are stored in the variable blade.materials. Each entry in
blade.materials is a MaterialDef object, which stores a name, elastic
properties, density, and strength properties among others (see :numref:`materialDefTable` in the :ref:`appendix`). It also stores the thickness that a single layer of that material in a composite is assumed to be, which can be important to know or edit when defining the thickness distribution of the blade’s components as just described.

After editing the design properties of a blade model as illustrated in
these few examples, a user should run the command

``blade.UpdateBlade()``

This function updates numerous internal private variables based on the
edited values in the public variables. Among other things, it
interpolates the properties that vary along the span of the blade to the
spanwise points specified in the variable ``blade.ispan``. These include
all the properties defined in ``blade.stations``, as well as the general
spanwise properties such as prebend, twist, etc. ``UpdateBlade`` also
updates the bill of materials for the blade, stored in ``blade.bom`` and
various details of the geometry, stored in ``blade.geometry``.

When the variables defining the blade design are set to satisfaction,
the blade object can be used to perform various operations for analysis
and optimization, such as generating representative structural models as
described in the next section.

.. _genBladeStructural:

Generating Representative Blade Structural Models
-------------------------------------------------

A NuMAD blade object can be used to construct structural models for
various types of analysis. Several tools exist that analyze
characteristics such as section stiffness, mass, and natural frequencies
of wind blades by representing them with low-fidelity beam models. These
include PreComp, BModes, and BPE. The most straightforward way of
invoking these capabilities is through the graphical user interface (for
details please see ref. [3]).

In addition to these, however, NuMAD 3.0 has many built-in functions for
performing high-fidelity analysis of a blade as a shell-element model in
ANSYS, which are easily invokable from a MATLAB script or command line.
These include analysis for maximum tip deflection, ultimate rupture
failure, global and local buckling, fatigue and natural frequencies and
are discussed in detail in :ref:`FEAOps`.

.. _AeroelasticSimRunIEC:

Aeroelastic Simulation and the ``runIEC`` Function
==================================================


A critical step in the design and optimization of any blade is
performing aeroelastic analysis to predict the behavior and response of
the blade under a range of expected wind and loading conditions. NuMAD’s
capability for performing this analysis is contained in the
``source\rotor_optimization\runIEC`` directory of the standard design codes
package. In the following sections the basic capability and
functionality of the ``runIEC`` package will be described, followed by
recent updates related to its operation from previous versions.

.. _useAndFunctOFrunIEC:

Use and Functionality of ``runIEC``
-----------------------------------

The main function called to perform aeroelastic analysis of a NuMAD
blade is

``source\rotor_optimization\runIEC\runIEC.m.``

This function calls on the accompanying tools in the directory to
process and output critical results, either by first running the
aeroelastic analysis for a set of requested simulation cases or by
processing a pre-existing set of output files from an analysis that was
previously run. The aeroelastic analysis is performed by calling FAST v7
software developed by the National Renewable Energy Laboratory
specifically for the analysis of wind turbines. FAST is a time-domain
solver, which output numerous variables defining the state of a wind
turbine at each time step for every simulation.

The call to ``runIEC`` takes three parameters, representing the
requested set of design load cases, a flag indicating whether to run the
aeroelastic analysis (as opposed to just processing an existing output
set) and a flag indicating whether to run in parallel, as shown below:

``output = runIEC(DLC,simflag,parflag)``

The inputs ``simflag`` and ``parflag`` should be set to 1 for yes and 0
for no. The design load case list ``DLC`` should be a cell array of
strings indicating the load cases desired. The supported load cases are
denoted by their codes in the IEC standard as described below:

**1.1:** Normal operating conditions with the turbine running and
connected to electric load. Simulations run with normal atmospheric
turbulence model, and 50-year-maximum loads are extrapolated based on
peak values in simulation results. Simulations are run for the range of
wind speeds and the number turbulence seeds specified in the variables
``params.ws`` and ``params.numSeeds`` in the ``runIEC_ipt.m`` file (``\examples\runIEC_ipt--EXAMPLE.m``).

**1.2:** Normal operating conditions with the turbine running and
connected to electric load. Simulations run with normal atmospheric
turbulence model, and fatigue damage is predicted based on cycle
counting of local peak values of loads in the simulation results.
Simulations are run for the range of wind speeds and the number
turbulence seeds specified in the variables params.ws and
params.numSeeds in the ``runIEC_ipt.m`` file. 

.. Note::
    DLC’s 1.1 and 1.2 call for the same simulation conditions, and if both are requested a single set of results is generated for both cases.

**1.3:** Normal operating conditions with the turbine running and
connected to electric load. Simulations run with extreme turbulence
model, and maximum loads taken directly from simulation results.
Simulations are run for the range of wind speeds and the number
turbulence seeds specified in the variables params.ws and
params.numSeeds in the ``runIEC_ipt.m`` file.

**1.4:** Transient modeling of an extreme coherent gust with wind
direction change with the turbine running and connected to electric
load, and maximum loads taken directly from simulation results.

**1.5:** Transient modeling of extreme wind shear with the turbine
running and connected to electric load, and maximum loads taken directly
from simulation results.

**5.1:** Modeling of emergency shutdown, maximum loads taken directly
from simulation results.

**6.1:** Modeling of extreme wind speed with the turbine in a parked or
idle state. Simulations run with extreme turbulence model, maximum loads
taken directly from simulation results.

**6.3:** Modeling of extreme wind speed and extreme yaw misalignment
with the turbine in a parked or idle state. Simulations run with extreme
turbulence model, maximum loads taken directly from simulation results.

As an example, the following call would run the FAST simulation for
cases 1.1, 1.3, and 6.1 without using parallel processing:

``output = runIEC({‘1.1’,‘1.3’,‘6.1’},1,0)``

The output of the ``runIEC`` function is a data structure reporting a
compilation of critical values from the results of all the simulations
run. These are mainly extreme values of structural quantities such as
reaction forces, reaction moments, deflections and strains at certain
points along the span of each blade. The specific simulation run and
time step at which these extreme values occur are also reported. In
addition, ``runIEC`` can also extract statistical information about
local peak values of these same structural quantities throughout all the
time histories for the purpose of predicting fatigue and long-term
extrapolated maximum values for DLC cases 1.1 and 1.2. The data in the
output structure as indicated in the example above is also written to a
spreadsheet file in the main blade model data directory titled
``IECDLC_Results.csv`` upon completion.

The runIEC function should be called from the main model data directory
of the blade to be analyzed. This directory must contain the necessary
FAST input files and airfoil data defining the model, a FAST output
directory and a NuMAD working directory. In addition, the model
directory must have a MATLAB script file with the name ``runIEC_ipt.m``,
defining a set of general simulation parameters referenced throughout
the process. A sample script showing the parameters that need to be
defined can be located at 

``\examples\runIEC_ipt--EXAMPLE.m``.

.. _recentUpdatesRunIEC:

Recent Updates to runIEC Functionality
--------------------------------------

Two main updates have been made recently affecting the ``runIEC``
process: 1) A change in the method of obtaining 50-year extrapolated
maximum values in the analysis of DLC case 1.1, and 2) An introduction
of tools for reading input files for OpenFAST. These changes are
described in the following two sections.

.. _methodsForDLC1p1:

Methodology for Analysis of Design Load Case 1.1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As per the IEC standards, an analysis of DLC case 1.1 involves the
simulation of a wind turbine under normal operating conditions for a
range of wind speeds sufficiently representing what is expected under
these conditions. Local peak values of important quantities such as
reaction loads on the turbine throughout the simulation time histories
are identified and used to derive statistical information. This
information is then used to predict 50-year extrapolated extreme values,
that is the highest/lowest value expected to be encountered over a
period of 50 years of normal operation, for each important quantity.

Needless to say, it is not common practice to actually simulate 50 years
of operation time for a turbine. The basic assumption is that if enough
time is simulated to extract reliable probability distributions of the
local peak values of the quantities of interest, those probabilities
would scale proportionally into significantly longer periods of time.
If, for example, an event is 1% likely to occur once over a period of 10
minutes, then it should be 10% percent likely to occur once over a
period of 100 minutes under the same conditions.

Following this reasoning, a load/quantity that has a probability of
3.805 × 10\ :sup:`-7` to exceed a certain value in a period of 10
minutes should be 100% likely to exceed that value once in 50 years of
operation. That value is the extrapolated 50-year extreme value, which
is the main result of interest for DLC 1.1.

While this conceptual approach is standard, there are several ways to
implement the 50-year extrapolation. The probability that the highest
local peak value :math:`F_{\text{ext}}` of a quantity exceeds a value
:math:`F` within a time period :math:`T` can be expressed mathematically
as follows:


.. math:: \text{Prob}\left( F_{\text{ext}} \geq F\  \middle| T \right) = P_{e}(F,T) = \int_{V_{\min}}^{V_{\max}}{\text{Prob}\left( F_{\text{ext}} \geq F \middle| V,T \right)p(V)\text{dV}}
    :label: maxPeak 

Where
:math:`\text{Prob}\left( F_{\text{ext}} \geq F \middle| V,T \right)` is
the probability of :math:`F_{\text{ext}}` exceeding :math:`F` within the
time period :math:`T` at a given wind speed :math:`V`, and :math:`p(V)`
is the probability density function of the wind speed. The
velocity-specific probability function is commonly computed as

.. math:: \text{Prob}\left( F_{\text{ext}} \geq F \middle| V,T \right) = 1 - \left( \text{CD}\left( F \middle| V,T \right) \right)^{N}
    :label: velProb

Where :math:`\text{CD}\left( F \middle| V,T \right)` is the cumulative
distribution function of the quantity for a given velocity :math:`V`
over the time period :math:`T`, and :math:`N` is the number of local
peak values of the quantity within the period :math:`T`. In words, Eq.
:eq:`velProb` is simply stating that the probability that the highest peak value
:math:`F_{\text{ext}}` exceeds :math:`F` is one minus the probability
that none do. This is a convenient way to evaluate the probability since
it is put in terms of the standard cumulative distribution function,
which is commonly available for most standard distribution types.

The task, then is to find the root :math:`F` of the equation

.. math:: P_{e}(F,T) = \frac{T}{T_{50\text{yr}}}


Or, for the typical 10-minute simulation time,

.. math:: P_{e}(F,T) = \frac{10}{60 \times 24 \times 365 \times 50} = 3.805 \times 10^{- 7}
    :label: prob10

and determine the 50-year extrapolated value for :math:`F`. One approach
to finding the root of Eq. :eq:`prob10` is to derive a single probability curve
representing :math:`P_{e}(F,T)` by compiling the local peaks from all
simulations throughout every wind speed, and fitting a set of parameters
defining a generalized extreme value distribution to the compiled data.
This curve can then be extrapolated beyond the range of data values to
obtain the root. This was the original approach used in the development
of the ``runIEC`` function. It proved to be potentially problematic,
however, in that the algorithm used to obtain the parameters for the
generalized extreme value curve often failed to converge reliably, and
when it did there were cases when it was infeasible to extrapolate the
curve to the range of the 50-year limit.

Consequently, it was determined that an alternative, more robust
implementation for finding the root of Eq. :eq:`prob10` was desirable, and new
procedure was developed and implemented for NuMAD v3. In the new
approach, a *set* of normal Gaussian probability distributions is
obtained, one for each wind speed for each quantity of interest. Those
distributions are used to evaluate the probability of exceedance (Eq.
:eq:`velProb`) for any arbitrary :math:`F` at a given velocity :math:`V`. Then
:math:`P_{e}(F,T)` (Eq. :eq:`maxPeak`) is evaluated with trapezoidal numerical
integration, using the appropriate distributions at each velocity, to
obtain the composite probability of exceedance of a quantity over
:math:`T` for an arbitrary :math:`F`. Finally, the root of Eq. :eq:`prob10` is
obtained using the bisection method, evaluating :math:`P_{e}(F,T)`
iteratively and converging to the solution for the 50-year extrapolated
:math:`F`. An algorithmic summary of the process is as follows:

For every quantity for which a 50-year extrapolated value is of
interest:

1) Process the aeroelastic simulation output files to extract the local
   peak values of the quantity at each simulated wind speed.

2) Calculate the probability distribution parameters (mean and standard
   deviation for normal Gaussian) to define the cumulative distribution
   function :math:`\text{CD}\left( F \middle| V,T \right)` for each wind
   speed.

3) Perform a bisection root-finding solve to find :math:`F` in Eq. :eq:`prob10`,
   each iteration evaluating :math:`P_{e}(F,T)` with trapezoidal
   integration as


.. math:: P_{e}(F,T) = \sum_{i = 1}^{N_{\text{ws}} - 1}{\frac{1}{2}\left( \left( 1 - \text{CD}\left( F \middle| V_{i},T \right)^{N} \right)p\left( V_{i} \right)\  + \ \left( 1 - \text{CD}\left( F \middle| V_{i + 1},T \right)^{N} \right)p\left( V_{i + 1} \right) \right)\Delta V} 


This is an improvement in robustness in the new approach, stemming from
two main aspects. First a normal Gaussian distribution is defined by the
mean and standard deviation, which can be directly and reliably computed
for any data set without any risk of convergence failure. Once those
parameters are known, the probability of exceedance can be extrapolated
to any value :math:`F`, without concern for the range of the original
data set. Second, the bisection algorithm for root finding is
fail-proof, provided that a suitable upper and lower bound is set in
which one root exists, and the function is continuous within those
bounds. It does not suffer from extreme gradient/slope values or
stability concerns.

.. Note:: 
    There can be slight differences between the 50-year extrapolation results obtained using different methods, and it is difficult to assert that any given approach is certainly more accurate or superior to another. The extrapolation process remains subject to further modifications and improvements moving forward.


.. _toolsForOpenFast:

Tools for Processing OpenFAST Input Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NuMAD is currently set up with FAST v7 as the default version of the
aeroelastic solver. Although FAST v7 is reliable and robust, it is now
an outdated version succeeded by the current equivalent, OpenFAST. The
framework of OpenFAST was designed to be customizable so that a broader
community of users could make contributions and variations for different
specific needs. As a result, the structure and organization of the input
files is significantly different from that of FAST v7. For the current
NuMAD release, some tools have been developed to process data from
OpenFAST input files to supplement what is currently there for FAST v7.

These tools primarily read different types OpenFAST input files, storing
the data in a MATLAB struct object which can then be edited and modified
for the purposes of design and optimization. Updated versions of the
input files can then be re-written from the modified data. The tools can
be found in the ``source\rotor_optimization\runIEC`` directory, along with
their FAST v7 counterparts.

As of the release of this document, OpenFAST remains in a state of
development, and new modules are coming online that will be increasingly
used in the future. The toolset for input/output processing is subject
to change to accommodate new input file formats, etc. moving forward.

.. _FEAOps:

Finite Element Analysis Operations
==================================

.. _meshGen:

Mesh Generation
---------------

The function called to generate the FE shell model in ANSYS of a NuMAD
blade is

``source\rotor_optimization\structOptimization\layupDesign_ANSYSmesh.m``

.. Note:: 
    It is currently necessary to have created the NuMAD input file from NuMAD 2.0, before attempting to run this function. See :ref:`NuMAD2p0` for further details. Thus, NuMAD 3.0 currently relies on NuMAD 2.0 to create an ANSYS mesh. As time and buget allow, the ``blade.generateFEA()`` will be fixed so that ``layupDesign_ANSYSmesh`` is not needed.

As an example, the following call would build a mesh for the blade blade
object

``layupDesign_ANSYSmesh(blade)``

This will generate an ANSYS models for NuMAD 2.0 input file named ``numad.nmd``.

If the NuMAD 2.0 input file is named differently then add a second argument as in:

``layupDesign_ANSYSmesh(blade,numadFile)``

where ``numadFile`` is a string with the desired filename. 


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
     The :math:`x_{i}` coordinate system is like the :math:`w_{i}` system but the first axis, :math:`w_{1}` points toward the leading edge instead of the flap direction.

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
    In the above equation, the minimum of :math:`M_{1}^{y}` is found instead of the maximum but :math:`\theta` still ranges from 0 to 180 deg. Unlike the resultants used for the deflection analysis which all occurred at a specific time in the OpenFAST simulations, each of the axial force, torsion, and bending moment resultants along the span could possibly come from different times. Thus, it is an artificial distribution suitable for design use. Also unlike the deflection case, the two bending-moment components (i.e. flap-wise and edge-wise bending) at a span location were projected onto 12 directions as defined by :math:`y_{1}` in :numref:`loadDirections`. All of the :math:`P_{i}\ `\ then transformed to the :math:`x_{i}` system in ``ad2ansys``.

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
is\ :math:`\ 1 \times n` MATLAB cell array where :math:`\text{\ n}` is
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
``ad2ansys.m``. ``ad2ansys.m`` generates a text file (usually called forces.src)
containing the APDL commands to apply nodal forces for every node on the
wetted area of the blade. Details of the approach are found in Ref. [1]
but with modifications to add axial loads. 

.. Note:: 
    The script assumes that the input forces are in the FAST coordinate system, :math:`w_{i}`, so a conversion is performed in the script to the ANSYS coordinate system, :math:`x_{i}`.

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
   :widths: 5, 1, 1
   :header-rows: 1
   

``meshFile`` is the name of the ANSYS model database to begin the analysis
from. It defaults to ``'master'``. analysisFileName is the name of the new
ANSYS model database which will store the analysis results. It defaults
to ``'FEmodel'``. np is the number of processors (CPUs) ANSYS should use for
each solve. It defaults to 1.

The output and it is a variable length struct. Depending on
which analysis flags are active, results can be accessed with

``result=layupDesign_ANSYSanalysis(blade,loadsTable,config)``

``result.globalBuckling``

``result.localBuckling``

``result.deflection``

``result.failure``

``result.fatigue``

``result.resultantVSspan``

``result.mass``

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
described as sandwich panels, local instabilities can be checked in
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

Make sure to set

``params.fatigueCriterion = 'Shifted Goodman'``

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
two fields

element1950: [1×1 struct]

element2558: [1×1 struct]

Looking at the first struct shows:

    >> myresult.element1950
         x3: [8×1 double] 
         eps11: [8×1 double] 
		 
         eps22: [8×1 double]
		 
         eps33: [8×1 double]
		 
         eps23: [8×1 double]
		 
         eps13: [8×1 double]
		 
         eps12: [8×1 double]
		 
         sig11: [8×1 double]
		 
         sig22: [8×1 double]
		 
         sig12: [8×1 double]
		 
         matNumber: [8×1 double]

where ``x3`` is the thickness coordinate, ``eps`` are the strains and sig are the
stresses. 

.. Note:: 
    Note the ``x3`` value that is zero corresponds to the reference surface in the shell model. 
	
In this example the shell offset in ANSYS was set to ``“BOT”`` so the zero appears the end:

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

    >> result.resultantVSspan{i}

where i is the i\ :sup:`th` load case in ``loadsTable``.

.. _optimization:

Optimization and Analysis Tools
===============================

The NuMAD package can be used to perform design optimization on a blade
structure. Through manipulation of the data defining the blade design
stored in the NuMAD blade object, a design can be tailored to best
satisfy a given objective, defined by a fitness or objective function. A
typical objective function may be the mass of the blade subject to
various constraints on aspects such as stress/failure, deflection,
buckling, fatigue or modal frequencies. This section contains some
general guidelines for setting up and running and optimization and
recent updates related to the process.

.. _runningOptimization:

Guidelines for Running Optimization
-----------------------------------

Optimization is simple in concept. In short it is the process of finding
the best possible state or performance of a system for what is invested
in it within allowable variable ranges and satisfying necessary
constraints. It is in no small way the essence of engineering and
design, but nothing could be more open-ended as there are virtually
limitless possibilities in how to define variables, objective functions
and constraints as well as what optimization algorithms settings and
parameters to use. A set of tools is offered for obtaining running
finite element, and other analysis to obtain the quantities needed to
evaluate the performance of a blade structure and set up an
optimization. Some general guidelines are given here to assist with the
process.

In addition, an example optimization script to demonstrate the application of the following concepts is included in 

``examples/optimization/exampleOptimizationDir``

This folder contains the main optimization script, ``optimizationExample.m``, the objective defined as a MATLAB function, ``objectiveExample.m``, along with a folder containing the ``.yaml`` file and loading data for the example blade.  The loading data is pre-generated using ``runIEC`` and the functions descriped in :ref:`AeroelasticSimRunIEC` and :ref:`FEAOps`.  An airfoil database directory is included as well, for reference in the NuMAD input file generation process.  To run the example script, place the exampleOptimizationDir folder with all its contents in a working directory of your choosing, and execute ``optimizationExample.m``.  The user is encouraged to read through the source code in main script and the objective function to understand the steps to the process, and the calls to NuMAD functions for various operations.  A good approach to putting together a customized optimization is to begin from these scripts and modify according to the specific needs at hand, while being mindful of the concepts presented in Sections :ref:`definingObjective` through :ref:`choosingOpimizationAlgor`. 

.. Note::
    Please update the ``ansysPath`` variable in the input settings to the actual path to your ANSYS executable.  


.. _definingObjective:

Defining the Objective Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first critical step in setting up an optimization is defining a
fitness, or objective function. This is a function that numerically
evaluates the favorability/suitability of a given design based on some
criteria. Any optimization process seeks to minimize or maximize its
objective function, and therefore it should represent whatever quantity
should be extremized to produce a favorable result. The objective
function must be defined as MATLAB function, which accepts a vector
representing the current values of the design variables, along with any
other inputs needed, calculates the value of the objective and returns
it as an output. For some optimization algorithms, the gradient of the
objective with respect to the design variables may be helpful to return
as well, though usually not strictly required. A common objective for
structural optimization would be mass, for example. If it is desired to
minimize the mass of a blade subject to some set of constraints, the
objective function should calculate the mass of the blade based on the
design variables and return it as the primary output. An example
objective function definition is given in ``examples/optimization/objectiveExample.m``, demonstrating
the use of NuMAD’s finite element analysis tools described in :ref:`FEAOps`.

.. _definingDesignVars:

Defining Design Variables
~~~~~~~~~~~~~~~~~~~~~~~~~

The next step is to define the design variables, or the variables that
can be changed in order to optimize the objective, and within what
ranges they can be changed. In NuMAD the design variables can, in
general, be any aspect of the blade’s design, as defined in :ref:`bladeDefAndTerms`,
in the description of the blade object. Typically, the design variables
should be taken from the input vector in the objective function and
assigned to their appropriate fields in the blade object, which can then
be used to evaluate any necessary quantities for the objective. The upper and lower bounds for each design variable are
usually passed to the main optimization function, which takes charge of
enforcing these bounds, but they may also be passed to the objective
function if it is useful. Determining appropriate initial values, upper
bounds and lower bounds can be challenging, and care should be taken to
make sure that the bounds do not violate physical limits, such as making
the spar cap wider than the blade or setting negative values for
composite layer thicknesses.

.. _definingConstraints:

Defining Constraints
~~~~~~~~~~~~~~~~~~~~

Next is to identify any constraints that must be imposed on the design.
A constraint is any condition that must be satisfied in order to make a
design feasible, or eligible as a solution. Constraints and objectives
are fundamentally similar, in that they are both essentially goals that
are hoped to be achieved in the final design. The primary difference is
that constraints typically impose a specific value, relative value or
range on a quantity which is not negotiable for a solution to be
acceptable. In contrast, an objective is a quantity that is sought to be
extremized, but the target value is not necessarily known, and it may be
compromised if necessary for the sake of constraints. There are
different types of constraints, and the best way to impose them can
depend on the type and the situation.

One type of constraint is a direct design variable constraint that
imposes a value or closed-form mathematical relationship between the
design variables themselves. The simplest example is the upper and lower
bounds placed to define the ranges of each variable. Usually these arise
from physical limitations and boundaries. For instance, if optimizing
the positions of two shear webs in a blade, the fore web must always be
positioned in front of the aft web, though bounds may still overlap.
This type of constraint be linear or nonlinear, and most optimizers can
easily accommodate them defined alongside of the objective function in
the form of a constraint matrix or similar structure, ensuring that they
are satisfied throughout the optimization.

Another type of constraint imposes a condition not on the design
variables themselves, but on some aspect of the performance of the
structure. An example of this is imposing that the maximum displacement
of the blade under an applied loading cannot exceed a certain limit.
Constraints like this are almost never linear or expressible in closed
form, and thus tend to be less straightforward to enforce.

One approach is to test the constraint at each design state throughout
the optimization, and simply discard any state that does not satisfy the
constraint, continuing to search in other directions of the design space
until a feasible state is found. This can be effective if the constraint
is not too restrictive, but it is potentially difficult, especially with
numerous complex constraints, to find states that fully satisfy all
constraints through blind trial-and-error. As a result, an optimizer can
run itself into the ground trying to find feasible candidates in this
manner.

An alternative method is to quantify the level of satisfaction of the
constraint in the form of a *penalty function*, and factor it in as part
of the total objective value. If the penalty function grows steadily
higher the further a constraint is violated, then the optimizer will
definitively favor solutions that satisfy the constraint and migrate
toward feasible solutions if it is violated. Returning to the former
example, if it is desired to ensure that the maximum displacement of a
structure under loading does not exceed a certain limit, a penalty
function could be derived in the form of, say, a power of the constraint
index as shown:

.. math:: U < U_{\max} \Longrightarrow \frac{U}{U_{\max}} < 1

.. math:: P = {c_{1}\left( \frac{U}{U_{\max}} \right)}^{c_{2}}


where :math:`c_{1}` and :math:`c_{2}` are some predetermined
coefficients. This method is generally effective, although it does not
guarantee that the final design perfectly satisfies the constraints, and
can require fine-tuning and adjustment of the constant parameters on the
part of the user. These are just some examples, and ultimately it is up
to the judgement of the user the most appropriate way to incorporate
their constraints in with the objective.

.. _choosingOpimizationAlgor:

Choosing an Optimization Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Numerous optimization algorithms exist, each with many options and
settings possible. A variety of functions are available built into
MATLAB, but many external packages are available as well. In choosing
the optimizer for a specific problem, it is helpful to classify
algorithms into two main types: gradient-based and gradient-free. These
two types have distinctly different strengths and weaknesses which are
important to understand in setting up a successful optimization.

Gradient-based optimizers work by beginning with the system and design
variables in a given initial state and taking a series of steps through
the design space searching for improved solutions. The search direction,
or the direction of the step by which the design variables are changed
each iteration is derived in some way from the gradient of the
objective, that is its derivatives with respect to the design variables.
The gradient is, by definition, the direction of most rapid
increase/decent for a function that is continuous and smooth throughout
the design space. The theory is that by using this direction, and
sometimes factoring in the Hessian of the objective and projecting out
the gradients of constraint equations, etc., it should be possible to
steadily reduce an objective from the initial state until it reaches a
minimum value.

Gradient-based optimizers are direct and efficient in suitable
applications. They typically require relatively few evaluations of the
objective function for an optimization, particularly if the
exact/analytical gradient is provided along with the objective. This
makes them attractive for problems with a large number of design
variables, or requiring high-fidelity analysis for evaluation of the
objective. They are, however, inherently local optimizers, meaning that
they seek out local extrema in the proximity of the initial state, which
may not be the most optimal solution in the whole design space. They
also require the objective function to be a continuous, smooth function
of the design variables which can be challenging to define in some
cases, while still representing the true objective that is sought after.
Examples of gradient-based optimizers include ``fmincon`` in the MATLAB
built-in optimization suite, or SNOPT, a widely-used sparse nonlinear
optimizer out of Stanford.

Gradient-free optimizers include a wide range of algorithms, which
search the design space in some way other than gradient-based. Because
of this, they do not rely on the objective having any particular
characteristics like continuity or smoothness, and they are in that
sense more robust and versatile than gradient-based alternatives. They
also tend to search the design space more thoroughly, although it should
be noted that it is rarely guaranteed that any optimization algorithm
will find the most optimal possible solution in the entire design space.
The methods that gradient-free optimizers use are diverse, but they
typically require a large number of objective function evaluations to
work effectively. If the objective is expensive to evaluate and the
number of design variables is very large, they can become prohibitively
costly, ineffective or both. In general, gradient-free methods are best
for problems with low to moderate fidelity in the objective simulations
and a few dozen design variables or less. Examples of gradient-free
optimizers include particle swarm and genetic algorithms (both offered
in the MATLAB suite), as well as some machine learning applications.

Keeping these basic guidelines in mind should assist in the process of
setting up an optimization for a specific application. In some
situations, it may be a prudent approach to first set up a gradient-free
optimization using lower fidelity analysis and a limited design space to
find a semi-optimal design, and proceed to fine tune it with a
higher-fidelity model and design space in a gradient-based optimizer.
Ultimately it is up to the user to determine the best course of action,
and the most important tool is your own judgement and ingenuity.


.. _recentOpiUpdates:

Recent Updates
--------------

Several updates have been made related to the tools for analysis and
optimization since the previous release of NuMAD. Among the most
significant is the addition of the finite-element-based analysis of
fatigue damage, and other module changes described in :ref:`FEAOps`. But in
addition to these, a few more subtle improvements have been made as will
be briefly described here.

For a period of time, when a finite element analysis of a blade was
performed based on the load output from a FAST/OpenFAST analysis, the
distributed load constructed and applied to the nodes of the model
consisted entirely of forces in the transverse (flap and edge)
directions, with no loading in the longitudinal direction. The reasoning
was that the flap and edge moments represented the most significant
loading on the blade for the purposes of predicting maximum stress, etc.
While this is generally true, it was decided that moving forward it
would be best to include the forces and moments in the longitudinal
direction to account for centrifugal effects and torsional moments for
the sake of completeness.

In the current version when forces are compiled from the FAST/OpenFAST
output in the functions:

``source\rotor_optimization\sim_tools\FastLoads4Ansys.m``

``source\rotor_optimization\sim_tools\getForceDistributionAtTime.m``

the longitudinal forces and torsional moments are compiled and applied
to the blade model along with the flap and edge moments. The appropriate
modifications were also made to the function:

``source\NuMAD_toolbox\ad2ansys.m``

to accommodate the longitudinal forces in the process. On a related
note, the forces and moments from the output files are given in a local
coordinate system at each point along the blade, which rotates along
with the structural twist defined for that point. For accuracy, the
forces and moments are now transformed to the blade global coordinate
system before being applied to the model. There is a new function
available to process a given output file, perform the transformation and
return the data in global coordinates:

``source\rotor_optimization\sim_tools\loadFASTOutDataGageRot.m``

Finally, when performing structural optimization, the blade model is
typically defined primarily by a ``.yaml`` file, which is read into an
instance of the NuMAD blade object for processing and design iteration.
But in several places throughout the process of ``runIEC``, and the
application of loads, information such as pre-bend, pre-sweep and
structural twist is taken from the FAST/OpenFAST files in the model
directory. To make sure the necessary information in these files is
consistent with that in the ``.yaml`` file, a convenient function was
built, named

``source\rotor_optimization\sim_tools\updateFASTFromBladeDef.m``

which updates and rewrites the fast files using the current data in a
blade file. This can be called immediately after reading the ``.yaml``
file, and before performing any analysis to ensure the consistency of
data.

.. _NuMAD2p0:

Coupling to NuMAD 2.0
=====================

Once the ``BladeDef`` object is defined, it is possible to interface
with prior versions of NuMAD which were GUI-centric. The function that
writes the input file to older versions of NuMAD from a blade object is

``source\preNuMAD\BladeDef_to_NuMADfile.m``

For example the function can be used as

``BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt')``

where ``blade`` is a blade object, ``'numad.nmd'``\ is the desired name
to be given to the NuMAD file and ``'MatDBsi.txt'`` is the desired
material data base name.


.. _GUI:

NuMAD GUI Mode
==============

In this version of NuMAD, the graphical user interface (GUI) can still
be accessed the same as it was in prior releases. Refer to the former user manual
(`SAND2012-7028 <https://energy.sandia.gov/wp-content/gallery/uploads/NuMAD_UserGuide_SAND2012-7028.pdf>`__) for
detailed instructions on how to use the GUI . For operating NuMAD exclusivley with the GUI refer
to the :ref:`intro-release-notes` on NuMAD v2.0.


.. _troubleshooting

Troubleshooting
===============

NuMAD will repeatedly print a message in the MATLAB command window to the effect of:

``Waiting for ANSYS to <`do something`>``
``Waiting for ANSYS to <`do something`>``
``Waiting for ANSYS to <`do something`>``

Common things to check are the ANSYS `.err` and `.log` for clues. 

.. _KnownIssues:

Known Issues
============
-  Fix ``blade.generateFEA()`` to create an ANSYS mesh without needed to use NuMAD 2.0 functions. This will make ``layupDesign_ANSYSmesh`` obsolete.


Trailing Edge Issues

-  There is a need to allow for the existence of a nonzero trailing edge
   thickness. The only way to achieve this in the current release is to
   specify that the airfoil in question is a flatback.

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

.. _PotentialFutureCapabilities:

Potential Future Capabilities
=============================

-  Allow for the option to make the blade entirely of solid elements

-  Incorporate adhesive modeling

-  Allow for FEA without the need of commercial FE licenses

-  Probabilistic flaw distributions

-  Incorporate Progressive damage

-  In addition to the ANSYS interface, add capability for a user to use
   other commercial FEA codes such as Abaqus and/or Nastran

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
   design tool (NuMAD V2. 0) for wind turbine blades: user’s guide."
   *Sandia National Laboratories Technical Report, SAND2012-7028*
   (2012). 

   
   
