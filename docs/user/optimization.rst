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

In addition, an example optimization script to demonstrate the application 
of the following concepts is included in ``examples/optimization/exampleOptimizationDir``

This folder contains the main optimization script, ``optimizationExample.m``, 
the objective defined as a MATLAB function, ``objectiveExample.m``,
along with a folder containing the ``.yaml`` file and loading data for 
the example blade.  The loading data is pre-generated using ``runIEC`` 
and the functions descriped in :ref:`AeroelasticSimRunIEC` and :ref:`FEAOps`.  
An airfoil database directory is included as well, for reference in the 
NuMAD input file generation process.  To run the example script, place 
the ``exampleOptimizationDir`` folder with all its contents in a working 
directory of your choosing, and execute ``optimizationExample.m``.  
The user is encouraged to read through the source code in main script 
and the objective function to understand the steps to the process, and 
the calls to NuMAD functions for various operations.  A good approach 
to putting together a customized optimization is to begin from these 
scripts and modify according to the specific needs at hand, while being
mindful of the concepts presented in Sections :ref:`definingObjective` 
through :ref:`choosingOpimizationAlgor`. 


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
the use of NuMAD's finite element analysis tools described in :ref:`FEAOps`.

.. _definingDesignVars:

Defining Design Variables
~~~~~~~~~~~~~~~~~~~~~~~~~

The next step is to define the design variables, or the variables that
can be changed in order to optimize the objective, and within what
ranges they can be changed. In NuMAD the design variables can, in
general, be any aspect of the blade's design, as defined in :ref:`bladeDefAndTerms`,
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
Examples of gradient-based optimizers include fmincon in the MATLAB
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
output in the functions: ``source\optimization\sim_tools\FastLoads4Ansys.m`` 
and ``source\optimization\sim_tools\getForceDistributionAtTime.m``

the longitudinal forces and torsional moments are compiled and applied
to the blade model along with the flap and edge moments. The appropriate
modifications were also made to the function: ``source\toolbox\beamForceToAnsysShell.m``
to accommodate the longitudinal forces in the process. On a related
note, the forces and moments from the output files are given in a local
coordinate system at each point along the blade, which rotates along
with the structural twist defined for that point. For accuracy, the
forces and moments are now transformed to the blade global coordinate
system before being applied to the model. There is a new function
available to process a given output file, perform the transformation and
return the data in global coordinates: 
``source\optimization\sim_tools\loadFASTOutDataGageRot.m``

Finally, when performing structural optimization, the blade model is
typically defined primarily by a ``.yaml`` file, which is read into an
instance of the NuMAD blade object for processing and design iteration.
But in several places throughout the process of ``runIEC``, and the
application of loads, information such as pre-bend, pre-sweep and
structural twist is taken from the FAST/OpenFAST files in the model
directory. To make sure the necessary information in these files is
consistent with that in the ``.yaml`` file, a convenient function was
built, named: ``source\optimization\sim_tools\updateFASTFromBladeDef.m``
which updates and rewrites the fast files using the current data in a
blade file. This can be called immediately after reading the ``.yaml``
file, and before performing any analysis to ensure the consistency of
data.
