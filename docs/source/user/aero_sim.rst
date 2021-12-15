.. _AeroelasticSimRunIEC:

Aeroelastic Simulation and the ``runIEC`` Function
==================================================

A critical step in the design and optimization of any blade is
performing aeroelastic analysis to predict the behavior and response of
the blade under a range of expected wind and loading conditions. NuMAD's
capability for performing this analysis is contained in the
``source\optimization\runIEC`` directory of the standard design codes
package. In the following sections the basic capability and
functionality of the ``runIEC`` package will be described, followed by
recent updates related to its operation from previous versions.

.. _useAndFunctOFrunIEC:

Use and Functionality of ``runIEC``
-----------------------------------

The main function called to perform aeroelastic analysis of a NuMAD
blade is ``source\optimization\runIEC\runIEC.m.``

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
set) and an IEC input structure in the form of an IECDef object, as shown below:

.. code-block:: matlabsession

    >> output = runIEC(DLC,simflag,IEC)

The ``simflag`` input should be set to 1 for yes and 0
for no.  The ``IEC`` input can be generated with the IECDef constructor,
 with a provided input file containing the necessary data.  The design load
 case list ``DLC`` should be a cell array of strings indicating the load 
 cases desired. The supported load cases are denoted by their codes in the
 IEC standard as described below:

**1.1:** Normal operating conditions with the turbine running and
connected to electric load. Simulations run with normal atmospheric
turbulence model, and 50-year-maximum loads are extrapolated based on
peak values in simulation results. Simulations are run for the range of
wind speeds and the number turbulence seeds specified in the variables
``windSpeeds`` and ``numSeeds`` in the IEC input file.

**1.2:** Normal operating conditions with the turbine running and
connected to electric load. Simulations run with normal atmospheric
turbulence model, and fatigue damage is predicted based on cycle
counting of local peak values of loads in the simulation results.
Simulations are run for the range of wind speeds and the number
turbulence seeds specified in the variables ``windSpeeds`` and 
``numSeeds`` in the IEC input file.

.. Note::
    DLC's 1.1 and 1.2 call for the same simulation conditions, and if 
    both are requested a single set of results is generated for both cases.

**1.3:** Normal operating conditions with the turbine running and
connected to electric load. Simulations run with extreme turbulence
model, and maximum loads taken directly from simulation results.
Simulations are run for the range of wind speeds and the number
turbulence seeds specified in the variables ``windSpeeds`` and 
``numSeeds`` in the IEC input file.

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
cases 1.1, 1.3, and 6.1:

.. code-block:: matlabsession

    >> output = runIEC({'1.1','1.3','6.1'},1,IEC)

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
directory and a NuMAD working directory.  For more information, see the
 runIEC example in the ``examples`` directory.

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
3.805 x 10\ :sup:`-7` to exceed a certain value in a period of 10
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
    There can be slight differences between the 50-year extrapolation 
    results obtained using different methods, and it is difficult to 
    assert that any given approach is certainly more accurate or superior
    to another. The extrapolation process remains subject to further 
    modifications and improvements moving forward.


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
be found in the ``source\optimization\runIEC`` directory, along with
their FAST v7 counterparts.

As of the release of this document, OpenFAST remains in a state of
development, and new modules are coming online that will be increasingly
used in the future. The toolset for input/output processing is subject
to change to accommodate new input file formats, etc. moving forward.
