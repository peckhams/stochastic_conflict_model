Stochastic Conflict Model
========

This stochastic model simulates the triggering of conflict in areas where "unrest" is high, and its subsequent spreading and/or resolution based on measures of local and nonlocal connectivity.  This model has been written as a Python package that can be installed in a conda environment. 

There is also a Jupyter notebook (see notebook folder) which has a graphical user interface (GUI), with sliders and text boxes, for running the model.  Instructions for installing the package are included in the notebook.  Values from the GUI are written to the model's configuration file.  Model output can also be visualized within the notebook with built-in graphics routines.  The model can also be run at a Python command line or at a Unix prompt.

This model has not yet been calibrated against conflict data (e.g. number of fatalities, type of conflict) because of the difficulty in obtaining reliable data that is grid-based vs. aggregated for admininstrative regions.  We are investigating various data sources including ACLED and GDELT for this purpose.

This model is partially compliant with the Basic Model Interface (BMI).

The main model parameters are as follows:

*conflict emergence factor*
This is a proportionality factor that is multiplied by the normalized unrest (U) to get the probability that conflict emerges in a given grid cell in the current timestep.  A real-valued number between 0 and 1.  Can be set with a slider in the GUI.
 
*conflict local spreading factor*
This is a proportionality factor that is multiplied by the normalized local connectivity (C1) to get the probability that conflict emerges in a given grid cell in the current timestep.  A real-valued number between 0 and 1.  Can be set with a slider in the GUI.

*conflict resolution probability*
The probability that the conflict in any grid cell will be resolved in the current timestep.  A real-valued number between 0 and 1.  Can be set with a slider in the GUI.

*unrest grid*
A spatial grid of values (2D array) that gives a measure of the "unrest", or potential for conflict to emerge, in each grid cell.  It depends on indicators such as population count, average rainfall rate, and many others.  This grid is currently pre-computed for the user as a function of several indicators.  It is fixed for each model run for a given region.

*local connectivity grid*
This is a spatial grid of values (2D array) that gives a measure of the "local connectivty", or potential for conflict to spread, in each grid cell.  It depends on indicators such as accessibility, road density, etc..  This grid is currently pre-computed for the user as a function of several indicators.  It is fixed for each model run for a given region.

*nonlocal connectivity grid*
This is a spatial grid of values (2D array) that gives a measure of the "nonlocal connectivty", or potential for conflict to spread to distant grid cells, in each grid cell.  It depends on indicators such as internet and cell phone access.  This grid is currently pre-computed for the user as a function of several indicators.  It is fixed for each model run for a given region.
