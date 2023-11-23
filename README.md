# Combinatorial-Vortex-Detection-Algorithm
A MATLAB Script for vortices in experimental velocity fields measured using particle image velocimetry (PIV) with a novel combinatorial vortex detection (CVD) algorithm. 

The CVD Algorithm implements an automated workflow for vortex identification and characterization in fluid flow data. It combines three techniques:

	1. Maximum Vorticity (MV) - Detects vortex candidates based on local vorticity maxima
	2. Cross-Sectional Lines (CSL) - Refines vortex core locations
	3. Winding Angle (WA) - Confirms and characterizes vortices using winding angle thresholding

This combinatorial approach enables robust vortex detection across varied flow conditions.

Usage:
The main CVD algorithm script sequentially calls the MV, CSL, and WA functions:

	1. Load input vector field data
	2. Set vortex detection parameters
	3. Initialize data storage
	4. Loop through vector field frames
		- Apply MV method to identify vortex candidates
		- Use CSL to refine vortex core locations
		- Perform WA analysis to confirm and characterize vortices
		- Store vortex data for each frame
	5. Generate visualizations and analyze results

Code Structure:
	- CombinatorialVortexDetection.m: Main CVD algorithm script
	- MaximumVorticityMethod.m: MV vortex candidate detection
	- CrossSectionalLinesMethod.m: CSL core localization
	- RegionsOfInterest.m: This function defines regions of interest (ROIs) around potential vortex cores
	- WindingAngleMethod.m: WA vortex confirmation/characterization
	- Circulation.m: Vortex circulation computation
	- Visualization functions: plotProfiles.m, plotStreamlines.m, plotStreamlinesWA.m, and plotVorticesMap.m

Requirements:
	- Image Processing Toolbox
	- Signal Processing Toolbox
	- Wavelet Toolbox
	- PIVMat Toolbox

References:
	
	- R.C. Strawn, D.N. Kenwright, J. Ahmad, Computer Visualization of Vortex 
	  Wake Systems, AIAA Journal. 37 (1999) 511–512. https://doi.org/10.2514/2.744.

	- H. Vollmers, Detection of vortices and quantitative evaluation of their 
	  main parameters from experimental velocity data, Meas. Sci. Technol. 12 (2001) 1199. 
	  https://doi.org/10.1088/0957-0233/12/8/329.

	- L.M. Portela, Identification and characterization of vortices in the 
	  turbulent boundary layer, Ph.D., Stanford Uni-versity, 1998.

Authors: Mathew Bussière, Guilherme Bessa, Bob Koch, and David Nobes.
         Department of Mechanical Engineering, University of Alberta, 
         Edmonton, Alberta

Contact: dnobes@ualberta.ca

Version: 1.0

Date: 10/6/2023
