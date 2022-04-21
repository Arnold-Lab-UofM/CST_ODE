This directory contains code and documentation of the Kirschner lab process for
performing sensitivity analysis of an ODE model using LHS (Latin Hypercube
Sampling) for generating input parameter data sets to efficiently sample the
model parameter space, run the model for each generated parameter set and PRCC
(Partial Rank Correlation Coefficient) to analyze the run results.

The LHS/PRCC method for performaing sensitivity analysis is from "A Methodology
For Performing Global Uncertainty And Sensitivity Analysis In Systems Biology",
Marino, et.  al., Journal of Theoretical Biology, 2008-09-07, doi:
10.1016/j.jtbi.2008.04.011. Another, prior, source is the book "Sensitivity
Analysis" by Saltelli, et. al., published by Wiley & Sons, ISBN: 0-471-99892-3.

The Kirschner lab uses Unix systems, Linux and Mac OS, including use of a
command line terminal window. If you use MS Windows a different approach may be
necessary. The documentation below on creating and running an ODE model assumes
a Unix like terminal window is available.

Version control is a necessity for any software process. If you are already
using version control, this directory should committed to your repository.

If you are not yet using version control I recommend using Git. This typically
comes preinstalled on most Linux and Mac systems (on Macs it may require
installing the Xcode command line utilities).

See git-scm.com for downloading git. See git-scm.com/book/en/v2 for an
excellent free online book on using git, including an overview of basic version
control concepts.

The documentation below does not use version control commands since your
version control setup probably will be different than the setup in the
Kirschner lab. It should be straightforward to adapt the procedures described
to an environment with version control.

************************************
Contents of Directory matlab-ode-lhs
************************************

Many of the files have "_new" in the file name. This is to distinguish this
version of Kirschner lab Matlab ODE LHS code from an older version. The
Kirschner lab keeps the older version in case older models need to be rerun, to
replicate prior runs whose results were used in publications. Only the new
version code is included here.

readme.txt
	This documentation file.

startup.m
	An example Matlab startup script for putting directory matlab-ode-lhs in
	the Matlab path.

lhs_ode_predator_prey
	A sub-directory with an example Lotka-Volterra predator prey model.

	This is used as a template for creating new models.

	It has 2 files:

		lhs_ode_predator_prey_ode.m
			This is a standard Matlab ODE model file that can be passed to a
			Matlab ODE solver.

			A copy of this file is copied and edited when creating a new model.

		lhs_ode_predator_prey_settings_new.m
			A settings file for running the model. The name of this file (or of
			a copy of this file, for a new model) is passed to script
			lhs_ode_run_new.m

			Here are the main items specified in the settings file. See a
			settings file, especially the comments at the start of the file
			for a complete description.

				The number of model runs to perform.

				The time points to save results for.

				Whether or not to perform PRCC analysis, and if so for which
				time points. These must be the time points to save results for
				or a subset of them.

				The name of the file containing the model equations, ex.
				lhs_ode_predator_prey_ode.m or a copy of that file for a new
				model.

				The parameter definitions including any parameter ranges to
				use for generating parameter samples by the LHS process.

				initial conditions, which may also contain LHS ranges.

				Optional pulse specifications. These may also contain ranges.

lhs_ode_run_new.m

	This is the main function to invoke in the Matlab command window to run a
	model. It takes one argument, the name of a settings file (without the
	".m"). Ex. lhs_ode_run_new('lhs_ode_predator_prey_settings_new')

check_function_handle.m
is_numeric_min.m
lhs_ode_default_output_labels_new.m
lhs_ode_define_run_matrices_new.m
lhs_ode_get_run_settings_new.m
lhs_ode_norm_new.m
lhs_ode_prcc_new.m
lhs_ode_unif_new.m
ltqnorm.m
	Support functions called, directly or indirectly, by lhs_ode_run_new.m.
	These should never need to be edited.

********************
Matlab starup.m file
********************

Directory matlab-ode-lhs needs to be in the Matlab path for this code to work,
since some of the functions call other functions. Rather than doing this
manually each time Matlab is run it is easiest to edit the startup.m file and
then move it the location where Matlab expects startup files to reside.

On Unix systems, such as Linux and Mac OS, this is
$HOME/Documents/MATLAB/startup.m. MS Windows may have a different arrangement
for this, for example a Windows registry entry.

************************
Creating a New ODE Model
************************

Copy the sub-directory for the example predator/prey model,
lhs_ode_predator_prey, to a new directory, with a descriptive name relevant for
the new model being created.

Rename the Matlab model and settings files in the new sub-directory to file
names relevant for your model.

For example, suppose we need th create an ODE model of TGFB and PGE2 receptor
dynamics of fibroblast and epithelial cells. "$" is the terminal window prompt
(which may be differnt on your system).

$ cd ~/models
$ mkdir tgfb-pge2-co-culture
$ cp -r matlab-ode-lhs/lhs_ode_predator_prey tgfb-pge2-co-culture
$ cd tgfb-pge2-co-culture
$ mv lhs_ode_predator_prey_ode.m lhs_ode_tgfb-pge2_ode.m
$ mv lhs_ode_predator_prey_settings_new.m lhs_ode_tgfb-pge2_settings_new.m

Edit the model file, lhs_ode_tgfb-pge2_ode.m, replacing the predator/prey
equations with the equations for your model and edit the settings file,
lhs_ode_tgfb-pge2_settings_new.m, replacing the predator/prey settings with
settings for your model - the parameters, initial conditions, etc.

See the comments at the start of the settings file, and at various other key
points in the file, for information on how to edit it.

Creating or updating the model and settings files is generally straightforward.
One common pitfall is specifying incorrect values for parameters and initial
conditions in a settings file. It is easy to make a typo when entering these
values. Also be sure the units are correct and consistent among the factors and
terms of a model equation. For example if one term of an equation is in
pico-grams/milli-liter and another term is in micro-grams/liter the results for
that equation won’t be correct.

***********************************
Running an ODE Model using LHS/PRCC
***********************************

After editing the settings file and/or model file for a model use script
lhs_ode_run_new to run the model.

It is best to have the Matlab current directory to be the directory with the
model file and settings file, and it is best to use a separate directory for
them, rather than put them in a directory with other unrelated files.

For example, in the Matlab command window
lhs_ode_run_new('lhs_ode_tgfb-pge2_settings_new'). Note that the ".m" is not
used when invoking the run script nor when specifying the settings file name.

The settings file contains a setting for the model file name, so there is no
need to specify that separately when invoking the run script.

The run script idoes the following

	* Read the settings file - actually run it as a Matlab function

	* Perform an extensive set of sanity checks on the settings. If any errors
		are encountered they are reported and the script exits.

	* Perform the LHS process on any ranges (for parameters, initial conditions
		or pulse specifications) to generate a matrix of parameter sets, with
		each row of the matrix a parameter set to use for a model run.

	* Runs the model for each generated parameter set, saving the results for
		that run in a data structure. The results are saved for the time points
		specified in the settings file.

	* If specified in the settings file performs a PRCC analysis on the results,
		for the PRCC time points specified on the settings file.

	* Prompts for optional saving of all the data structures produced.
