Merge TimeSteps Plugin
======================

Add the filter 'Append Attributes Over Timesteps' (vtkMergeArraysAndTimeSteps)
that is based on the 'Append Attributes' (vtkMergeArrays) filter. Arrays from each
timestep from each input are append to the output, with the name
'originalName_input_#N_ts_#TS' where #N is the input index and #TS is the
timestep.
