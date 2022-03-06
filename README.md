# Single Photon Timing Analysis and Low-Gain Photon Timing Calibration Code
Developed by: Alan Kahn  (alan.kahn@cern.ch), ~2018

This code needs ROOT to be setup in order to run:

```
setupATLAS
lsetup root
```

This directory contains two code packages:

1. Timing Analysis
2. Calibration Analysis

## Timing Analysis

This code is to make analysis plots similar plots to those found in the Run 1 paper, i.e. Delta Z, t_gamma, etc.

- Run/NtupleCode/TimingAnalysis.C

Run with ``` ./runTimingAnalysis.sh ``` in the Run/ directory.

Output is saved to Run/NtuplePlots/exot6_plots.root

## Calibration Analysis

This code is to calibrate the timing of low gain photons in the 2016 data. 

- CalibrationAnalysis/Root/CalibrationAnalysis.C

Run with ``` ./runCalibrationAnalysis.sh ``` in the Run/ directory.

Output is saved to Run/NtuplePlots/calibration_exot6.root

If we plot the time of low-gain photons right out of the box, we get a distribution that is off-center and wide.
	Ideally, we want a distribution that is both centered at 0 and narrow. In order to do this, we have to take our
	original distribution, plot it vs a variable that we think could correlate with the time value in some way, and
	determine a correction so that once applied, the distribution becomes as close to a flat line at 0 in this 2-D
	plot. Then, we look at the 1-D corrected distribution, see how it was affected, and repeat the process again with 
	more and more variables until we're satisfied with the combined result from all corrections. 

Thus, this procedure is done in passes, with each pass providing a correction. This also means that this cannot be
	done all at once; rather, we need to go through all events for the first pass, determing a correction, modify the
	time by this corrected amount, and then go AGAIN through all of the events, now using their corrected time from 
	the previous pass, to determine another correction, etc.

In the code, you will see multiple event loops with a comment just before them stating which pass the loop corresponds to.

- **Pass 0**: This applies a time-of-flight correction. When an ECAL cell registers a hit, the time it outputs is calculated as if an object
	came from the center (z=0 point) of the detector. However, the primary vertex is not always at the very center; rather it is distributed
	narrowly around it. The time-of-flight correction takes into account the difference in path length, and therefore time, that arises
	due to a photon coming from the primary vertex rather than the center. I call this pass 0 since there is no correction that 
	needs to be found, it is just a correction from physics that needs to be applied to all objects.

- **Pass 1**: This correction is determined from a profile plot of pass 0 time vs run number. The profile plot is 2-dimensional, and is binned
	along the X-Axis, where each bin corresponds to a particular run number. The Y-Axis is the mean time of photons that are in the X bin 
	that it corresponds to. If we make this plot, we see that each run number corresponds to a different average time. Thus, each Y-value is
	the amount that all the photons' times should be subtracted by to arrive at a distribution centered around zero - i.e., the means for each
	run number are shifted to be zero. The before and after plots for this correction can be seen by opening up the root file in 
	Run/NtuplePlots/calibration_exot6.root and looking at the p_t_pass0_vs_runN and p_t_pass1_vs_runN plots. The corresponding 1-D distributions
	are h_t_pass0 and h_t_pass1.
	NOTE: The ECAL is separated into 4 partitions. 2 comprising the barrel, and 2 comprising the endcaps. They are split in the middle of the detector
	and are called EMB A and C, and EMEC A and C. In order to be more accurate, we make this correction for each partition. So in total we have 4
	profile plots that we determine corrections for, 1 for each partition.

- **Pass 2**: This correction is just like pass 1, however, instead of plotting time vs run number, we plot it vs front-end-board (FEB) number.
	Only one correction has to be made for each FEB, rather than splitting it up into partitions like pass 1, so there's just one plot as a result.
	You can see the before and after with the p_t_pass1_vs_febN and p_t_pass2_vs_febN plots. Corresponding 1-D plots are h_t_pass1 and h_t_pass2.

- **Pass 3**: If you look at the p_t_pass_2_vs_febN plot, you'll see that the first half of it looks like alternating narrow and wide distributions.
	The width of these regions in X is 32, meaning that each 32 FEBs behaves differently. As a result, we define a FEB SLOT as a region of 32 febs.
	The first 8 slots correspond to the barrel. In order to correct for this strange behavior, we decide to plot the pass 2 time vs the ECAL cell
	energy for each slot. These plots are listed as p_t_pass2_vs_E_cell_slot_# (# = 0->7). We see that each distribution has some systematic
	curve. We then fit this curve to a polynomial, determine the equation of fit, and use this equation to correct our time again. This is only 
	done for the first 8 slots, which means that this correction is only applied to photons in the barrel, and not the endcap. The corrected profile
	plots are p_t_pass3_vs_E_cell_slot_# and the respective 1D plots are h_t_pass2_slot_# and h_t_pass3_slot_#.

- **Pass 4**: This is like pass 3, in the sense that the correction is determined by an equation of fit. The variables we are plotting against here are
	f1 and f3, the fraction of energy deposited in the first and third layer of the ECAL cell respectively.

- **Validation Loop**: This applies all corrections on a different set of photons than the corrections were determined for. It's here to see how well
	these corrections work on an independent set of photons to see if we get comparable results to our control sample, the one we determine the 
	corrections from. Events in the control region and validation region are differentiated by their values of MET for the whole event. Events
	with MET < 100GeV correspond to the conrol region, while events with 100GeV < MET < 200GeV correspond to the validation region. Events
	with MET > 200GeV correspond to our signal region, which we will not look at until we have our analysis all ready and approved for unblinding.

## Plotting Macros

I also included some plotting macros in the Run/NtupleCode directory, listed as makePNGs.C and fit.C

makePNGs.C opens up a file, chooses a histogram or a profile plot, sets some of the plotting parameters to make it look nice, and then saves
	it to at PNG in the Run/NtupleCode/Plots folder. Feel free to copy any of the code snippets and change it to match what you'd like to plot.

fit.C is what I use to plot equations on profile plots. A tutorial of fitting in root can be found here https://root.cern.ch/root/HowtoFit.html
	Output from this code is also saved in Run/NtupleCode/Plots

