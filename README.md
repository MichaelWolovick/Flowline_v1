# Flowline_v1
This is the code for the flowband model used in Wolovick and Moore, 2018

These are the matlab scripts for the flowband model I used to study glacial geoengineering in this paper:  https://doi.org/10.5194/tc-2018-95.  

The main script for the model is Flowline_v1.  All scripts are extensively commented, and that one should tell you how to run it.  In addition, FlowlineBundler_v1 can be used to run many models in a sweep through parameter space.  

The behavior of Flowline_v1 is controlled by two sets of inputs:  the parameters that are set in the script itself (in a section conveniently labelled "Parameters") and the parameters from the input file that define the initial conditions and boundary conditions.  In some cases, a boundary condition can be defined in either place; the switching between them is controlled by setting the value of the parameter in the script itself to the string 'file'.  Usually that is used in cases where I wanted to use a spatially or temporally varying BC instead of a constant value.  In all cases in which that is possible, the comment next to the parameter will say that that is an option.

In the case that FlowlineBundler is being used to run many instances of the model, then Flowline_v1 should be switched to function mode rather than script mode by uncommenting the function definition line at the top of the code.  Be aware that any variables that are imported from the top when in function mode need to be commented out where they are defined lower down in the script (matlab's editor should underline them indicating that their value is being overridden).  

Flowline_v1 is also capable of doing a number of other things that were not included in that paper.  It is based on a model that I wrote during my thesis for solving ice flow and deformation problems associated with basal melt/freeze and travelling slippery/sticky patches.  As such, this model can also compute ice temperature, it can use both higher order (Blatter-Pattyn) and Full Stokes velocity solvers, it can run a basal hydrology model, and it can advect tracers and track layers.  I've included the subscripts for that functionality here as well.  All of those features can be turned on or turned off using variables under the "Model Parameters" section.

I have included more scripts than are necessary here.  I basically just uploaded most of the scripts that I wrote while working on this project, so there are multiple copies of some scripts.  The version that is called from Flowline_v1 is generally the "official" version.

I have also uploaded the input files, produced by inverting the width-averaged velocity data.  Those are the .mat files.  The original input files included a whole ensemble population from the last generation of the evolutionary algorithm I ran to do the inversion, as well as all sorts of info about the inversion, but those files were too big to upload to github so I created trimmed down files containing just the "*_input" variables that are actually information needed to start the model.  The "_slim" suffix indicates that the files are the trimmed-down versions.  You should be able to start the model using those input files.  Flowline_v1 only loads variables matching the pattern "*_input" from the input files.

Some of these scripts are not part of the model per se, but I used them to analyze model output or to generate figures or animations for the paper.

Questions, comments, or requests for clarification can be sent to wolovick@princeton.edu or Michael.Wolovick@gmail.com.

Hope this is useful,

Mike Wolovick

ps.  Yes I know that this is actually a flowBAND model (it can do variable width) rather than a flowLINE model.  Sue me.
