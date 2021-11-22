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

Eample 2 Name
-----------------


Eample 3 Name
----------------- 