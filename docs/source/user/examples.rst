.. _examples:

Examples
========


This sections shows how to perform basic operations in a step-by-step manner. Refer to ``$NuMAD/examples`` for files needed for running each example. 

.. _excelToNumadExamples:
   
Excel To NuMAD Examples
----------------------------------------  

NuMAD has the capability to generate blade models from Excel-based input files. This capability is offered for both the GUI and blade object 
modes of operation. While both the GUI and blade object offer excel compatability, the excel input file read by the GUI is not compatable with 
the blade object nor is the blade object excel input file compatable with the GUI at this time. To ease use of the excel capabilities within 
NuMAD 3.0 examples for generating a blade model from excel for the GUI and the blade object are given. It is recommended if the user is going to 
make their own blades using the Excel input-file capability, to use the example files as a base and modify to minimize errors in formatting.
  
.. _excelToObject:

Excel To Object
~~~~~~~~~~~~~~~

The example files needed to begin generating blade models from excel in the blade object are located in folder ``examples\ExcelToObject``. In the folder
there is an airfoils folder with the necessary airfoils for the example blade, the excel input file ``Excel2ObjectExample.xlsx``, and the 
example driver script ``excel2NumadObject.m``. Running the driver script will generate a blade object from the excel file, produce an equivalent
beam model using NuMAD objects built-in ability to call PreComp, and finally generate an ANSYS shell model. 

.. _excelToGUI:

Excel To GUI
~~~~~~~~~~~~~~~

The example files needed to begin generating blade models from excel in the NuMAD GUI are located in folder ``examples\ExcelToGUI``. In the folder
there is the excel input file ``Excel2NumadGUIExample.xlsx`` and a bill of materials template (for tracking purposes) ``BOM_template.xlsx``. 
To begin generating blade models in the NuMAD GUI from an excel file, first open up the NuMAD GUI with the command ``numad``. Next, select 
File -> New Model from the drop down menu. Then, select File -> Save model, this will generate the necessary generic airfoil folder and ``MatDBsi.txt`` 
files needed. Finally, select File -> XLS-to-NuMAD, and select the excel input file in your case directory.

.. _runIEC:

Run Aeroelastic Simulations to Analyze IEC-Standard Design Load Cases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example script is included to demonstrate the use of the runIEC function, which performs aeroelastic simulations using FAST in order to analyze the performance of a blade design under the IEC-standard design load cases.  The script can be found in the directory ``examples\runIEC\IECLoadsAnalysis.m``.  The basic steps performed by the script are as follows:

    -Set basic user options for the run, including which design load cases to run, whether to run a FAST simulation or use existing output files, and whether to generate loads tables and run an optional ANSYS FEA analysis from the results.
	-Include the directory containing the NuMAD source code files in the MATLAB operating environment, and establish the paths for external analysis packages.
	-Load the data from the IEC input file into an IECdef class object.
	-Load the pre-saved blade object for the example blade into the MATLAB workspace.
	-Create a NuMAD input file from the blade object.
	-Execute the runIEC function based on the user options, storing the output structure in IECOutput
	-If the option is selected to run an ANSYS analysis:
	    -Create a working ANSYS subdirectory, and copy in the NuMAD input file and the airfoils directory.
		-Generate loads tables from the IECOutput data.
		-Set configuration variables for the ansys deflection, rupture/buckling, and frequency analysis
		-Execute ANSYS analysis and store output.
		
The directory ``examples\runIEC`` also contains various files and folders used in the analysis.  The 'AeroData' folder contains aerodynamic data used in the FAST analysis by AeroDyn.  The 'airfoils' folder contains airfoil shape profiles for the design.  The 'init' folder contains sample input files that are used as templates for the input file creation tools within runIEC.  The 'out' folder is the directory where FAST output files are to be stored.

Alongside the main IECLoadsAnalysis script are an example set of FAST input files defining the turbine structure and simulation parameters.  In addition, an input file defining settings for the runIEC function is included, with the example blade stored as a NuMAD blade object, and a file specifing information for the turbine control system, named ``pitch.ipt``.  The IEC input file must be provided by the user any time runIEC is executed, but the blade object is really only needed if it is desired to run a subsequent ANSYS analysis.  The user may find it useful to build cased by modifying the given example to become comfortable with the workflow.

Note that to execute runIEC in its full capacity, the user must have FAST, Crunch, TurbSim and IECWind installed as additional packages, and the executable paths must be defined in the addNumadPaths script in the main NuMAD directory.  It can happen on occasion that a FAST simulation may terminate prematurely, leaving an output file without the entire requested time history of results.  This is usually due to a numerical aspect such as the turbulence seed that only applies to a given run, and it can cause the execution of Crunch to abort, and the runIEC to fail to run to completion.  Incidents such as this can be dealt with on a case-by-case basis, by either changing the seeds for the turbulent cases or adjusting the range of windspeeds in the IEC input file.  

The main results of runIEC are stored as IECOutput in the form of a MATLAB structure in the provided example, and also written in spreadsheet format as .csv files.  Thes results, and the subsequently generated loads tables can be used to assess performance metrics such as maximum tip deflection and failure for a blade, and can be used in optimizations, as detailed in the following example.

.. _Optimization:

Run Optimization to Refine the Design of a Blade to a Given Objective
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

An example script is included which demonstrates the basic process of setting up and running a blade structural design optimization problem within the NuMAD framework.  It is located in ``examples\optimization\optimizationExample.m``, alongside all resources needed to run it.  Optimization in general is a very open-ended process and the provided example is by no means a comprehensive demonstration containing every functionality or method that could be needed.  It should, however, provide a solid starting point and guide for integrating the analysis functions within NuMAD into an optimization process.  The main script walks through the following steps:

    -Set the basic user settings for the optimization, including directory of starting model data, and parameters for how to run the optimization.
	-Set configuration variables for the ansys deflection, rupture, buckling, fatigue and frequency analysis.
	-Include the directory containing the NuMAD source code files in the MATLAB operating environment, and establish the paths for external analysis packages.
	-Set up parallel working directories as appropriate.
	-Read the .yaml file defining the blade design and geometry into a NuMAD blade object, and read the IEC input file.
	-Load the pre-saved tables of applied loads into the MATLAB workspace.
	-Execute the optimization based on the method defined in user settings.  Included options are gradient based optimization (OptAl='gradient'), particle swarm (OptAl='particleswarm') or evaluate objective with no optimization (OptAl='objective').
	
Note that included are only two options for optimization algorithms, but MATLAB has many more which can be invoked similarly to those in the example.  See MATLAB documentation for more details.  For any optimization, the user must define an objective function which calculates and returns the value to be optimized.  In the provided example the objective function is the mass of the blade, with penalty constraints on maximum deflection, maximum material failure index, minumum buckling load factor, maximum fatigue damage, and natural flap frequency.  The design variables are defined to be the thicknesses of the individual blade components along the span.  The objective function is defined in ``examples\optimization\objectiveExample.m``, and walks through the following steps:

    -Determine the current parallel working task and change to the appropriate directory.
	-Modify the thicknesses of the blade components based on the values in the input design variable vector.
	-Create an ANSYS shell model of the blade.
	-Calulate the penalties due to the constraints using NuMAD's ANSYS analysis functions, storing the sum in the objective variable.
	-Add the total mass of the blade to the objective value.
	-Write results and key quantities to objective history file, and return the objective value to the optimizer.
	
The example script has 2 companion folders alongside it.  The ``airfoils`` folder contains a collection of airfoil profiles used by the model.  The ``exampleBlade`` folder contains data defining the given example model, including a .yaml file with all the blade geometry and material properties, MATLAB data files defining the applied loads to be considered, an IEC input file and a MATLAB data file containing the rain-cycle-counting data pertaining to fatigue analysis.  These can be generated from the runIEC function as demonstrated in the previous example.  After executing an optimization, if it runs to full completion, the ``examples\optimization`` directory will have a file containing a comprehensive list of the outputs at every design state encountered in the optimization process, named ``completeObjectiveHistory.txt``, as well as a file containing the final optimized values of the design variables, named ``FinalDVar.txt``.  The user has free range to modify this example as they will, to suit their individual needs.
	