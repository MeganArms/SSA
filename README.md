# SSA
**Single Spot Analysis**
Developed in The Hess Lab for Nanotechnology and Synthetic Biology by Megan Armstrong
Department of Biomedical Engineering
The Fu Foundation School of Engineering and Applied Science
Columbia University

## What is the purpose of this directory?
The files here are designed for the organization of data output from single molecule spectroscopy experiments analyzed with the impage processing MATLAB program developed by Raghuveer Parthasarathy (Nature Methods volume 9, pages 724–726 (2012)
doi:10.1038/nmeth.2071). ADD BG ON WHY AND WHY THE SF ISSUE 

Input: output from the RP analysis program.

Required: MATLAB R2018a and RP analysis program/GUI, available [here](http://pages.uoregon.edu/raghu/particle_tracking.html).

Output: Individual video files with variables **C**, **B**, **F**, **M**, **S**, **R**, **Frame**, **spotEvents**, **eventfreq**. Cell arrays with aggregated data otuput **survFuncs**, **intensities**, **numbers**

Definitions:

- **C**: cell array with each cell containing the `i,j` coordinates of the position of the particles on subsequent frames in row 1 and row 2, respectively. 
- **B**: cell array with each cell containing the integrated intensity of the particle on each frame that it appears.
- **F**: cell array with each cell containing the frame number on which the particle appears.
- **M**: cell array with each cell containing the particle identification number corresponding to the column in which the particle's parameters are in the objs_link output file from the RP analysis program.
- **S**: cell array with each cell containing the "sigma", or width, of the intensity spot indicating the presence of the particle.
- **R**: structure with one field where each row containing the particle identification number corresponding to the column in which the particle's parameters are in the objs_link output file from the RP analysis program. This may be eliminated in future versions. 
- **Frame**: cell array with each cell containing the particle identification numbers on that frame.
- **spotEvents**: structure with one field where each row containing the particle identification number after the trajectories have been connected across dark frames.
- **eventfreq**: Nx1 double with each row containing the number of events that occur on a single location of interest. N is the same length as `spotEvents`.
- **survFuncs**: cell arrray with each row containing collected information on survival functions from a video or set of two videos. Each column contains:
	* 1) survival function connected across dark frames, using the Kastantin-Schwartz correction method.
	* 2) time points for the above survival function
	* 3) standard error of the mean at each time point.
	* 4) sum of the correction factors for the above survival function
	* 5) numerator at each point for the above survival function, i.e. the corrected counts for each time point in the survival function.
	* 6) survival function connected across dark frames, using the Kaplan-Meier correction method.
	* 7) time points for the above survival function
	* 8) double array of residence times used to generate the survival function
	* 9) logical array indicating if the residence time was censored (1) or not (0)
	* 10) survival function using the Kastantin-Schwartz correction method.
	* 11) time points for the above survival function
	* 12) sum of the correction factors for the above survival function
	* 13) numerator at each point for the above survival function, i.e. the corrected counts for each time point in the survival function.
	* 14) survival function using the Kaplan-Meier correction method.
	* 15) time points for the above survival function
- **intensities**: cell arrray with each row containing collected information on intensities from a video or set of two videos. Each column contains:
	* 1) double array of the total integrated counts of each, unconnected trajectory
	* 2) double array of the number of frames upon which each unconnected trajectory appears
	* 3) double array of the number of bright events per trajectory, repeated for number of bright events in each connected trajectory.
	* 4) double array of the number of frames upon which each connected trajectory appears, repeated for the number of bright events in every trajectory.
	* 5) double array with the length of all of the dark events in every connected trajectory.
	* 6) double array with the average brightness of each connected trajectory
	* 7) double array with the range of brightnesses of each connected trajectory
	* 8) double array with the maximum brightness of each connected trajectory
	* 9) double array reverse cumulative sum of the average brightness on each frame over all trajectorie.
- **numbers**: cell arrray with each row containing collected information on particle counts from a video or set of two videos. Each column contains:
	* 1) double array with the first frame number for each conencted trajectory
	* 2) double array with the last frame number for each connected trajectory
	* 3) double with the number of trajectories appearing on the first frame
	* 4) double with the fraction of trajectories appearing on the first frame
	* 5) double with the number of trajectories appearing on the last frame
	* 6) double with the fraction of trajectories appearing on the last frame
	* 7) double with the number of connected trajectories
	* 8) double with the number of trjaecoties
	* 9) double with the number of objects detected
	* 10) Nx1 double with each row containing the number of events that occur on a single location of interest. N is the same length as `spotEvents`, or the sum of the length of `spotEvents` from two videos.
	* 11) 1x1 double with the number of frames in the video.

## Quick Start Guide to SSA

### Follow instructions in TrackingGUI_rp_Manual.pdf

### Load output file from RP analysis program 
- Either double click on the file, drag and drop the file into the command window, or run the following in the command window with the appropriate path and filename:
```
load('/Path/Filename.mat')
```

### Clean up the data
- Run
```
objs_link1 = objs_link;
clear objs_link objs
[C1,B1,F1,M1,S1,R1,Frame1] = collate(objs_link1);
[spotEvents1, eventfreq1] = eventlinks9(2,C1,M1,B1,F1);
save('/Path/Filename.mat')
clear
```
- If there are two videos for the same condition, load the second file and run
```
objs_link2 = objs_link;
clear objs_link objs
[C2,B2,F2,M2,S2,R2,Frame2] = collate(objs_link2);
[spotEvents2, eventfreq2] = eventlinks9(2,C2,M2,B2,F2);
save('/Path/Filename2.mat')
load('/Path/Filename.mat')
```

### Aggregate output for further analysis of two videos per condition
- Replace `LAG` with the time between frames, and `MMDD ExperimentalConditionX` with the date and experimental condition and run
```
[survFuncs,intensities,numbers] = generateoutput(LAG,Nframes,C1,C2,F1,F2,spotEvents1,spotEvents2,objs_link1,objs_link2,eventfreq1,eventfreq2);
save('MMDD ExperimentalConditionX output.mat','survFuncs','intensities','numbers')
clear
```

- If there are multiple experimental conditions to be compared, repeat the data clean up step and run
```
load('MMDD ExperimentalConditionX output.mat')
[survFuncs,intensities,numbers] = generateoutput(LAG,Nframes,C1,C2,F1,F2,spotEvents1,spotEvents2,objs_link1,objs_link2,eventfreq1,eventfreq2,survFuncs,intensities,numbers);
save('MMDD ExperimentalConditionX output.mat','survFuncs','intensities','numbers')
clear
```

### Aggregate output for further analysis of photobleaching experiments
- Replace `LAG` with the time between frames, and `MMDD ExperimentalConditionX` with the date and experimental condition and run
```
[survFuncsPB,intensitiesPB,numbersPB] = generateoutputPB(LAG,Nframes,C1,F1,spotEvents1,objs_link1,eventfreq1);
save('MMDD ExperimentalConditionX output.mat','survFuncsPB','intensitiesPB','numbersPB')
clear
```

- If there are multiple experimental conditions to be compared, repeat the data clean up step and run
```
load('MMDD ExperimentalConditionX output.mat')
[survFuncsPB,intensitiesPB,numbersPB] = generateoutput(LAG,Nframes,C1,F1,spotEvents1,objs_link1,eventfreq1,)
survFuncsPB,intensitiesPB,numbersPB);
save('MMDD ExperimentalConditionX output.mat','survFuncsPB','intensitiesPB','numbersPB')
```

### Photobleaching analysis
- Find the exponential decay of the photobleaching residence time curve and run `sfcorr.m` to remove the photobleaching effect.
```
load('MMDD ExperimentalConditionX output.mat')
mu = expfit(survFuncsPB{1,8},0.05,survFuncsPB{1,9});
sc = cell(size(survFuncs,1),1); s = sc;
for k = 1:size(survFuncs,1)
	sc{k} = sfcorr(survFuncs{k,7},survFuncs{k,6},mu);
	s{k} = sfcorr(survFuncs{k,15},survFuncs{k,14},mu);
end
save('MMDD ExperimentalConditionX output.mat')
```


## Potential output

### Survival function of connected trajectories
- Run
```
load('MMDD ExperimentalConditionX output.mat')
plotallcurves(survFuncs(:,7),sc(:))
```


