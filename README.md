# optCorr_script
  optCorr is a simulation & optimization based algorithm for natural isotope abudance correction. 
The state-of-the-art resolution-dependent algorithms like Accucorr and IsoCorv2 determine whether a non-tracer natural isotopic component should be corrected or not based on the MS resolution provided and define the proper correction matrix accordingly to correct only for those "unresolved" ones.  It is assumed that, the observed pre-correction signal is the simple summation of the unresolved peak abundences. However, this assumption may not be valid since peak top (instead of peak area) is often used as the signal readout when the profile spectrum is turned into centroid.  In fact, the unresolved contaminant peaks may only partially contributes to the observed signal, but are 100% subtracted in Accurcorr and IsoCorv2, which can cause over correction. optCorr uses a built-in MS-iso simulator that simulates the MS profile spectrum of an arbitrary labeling% at a given resolution, from which a centroid spectrum and isotopic distribution signals are extracted and compared to the experimental. the right amount of correction is determined by simulation. True labeling% is obtained by non-linear optimiation that minimized the least squre error between the simulated and observed signals.

optCorr currently supports single labeled 13C, 15N, 2H or 18O isotope correction. 

<br />  This is a light script verion of optCorr that contains the core "optcorr" function and a wrapper for batch processing elmaven output files:
<br />   function [distout,err,ysim]=optcorr(mz,atoms,type,abd,imp,fwhm,dmz,distin,option)
<br />   function optcorr_batch(varargin)  
## Requirement
 <br /> Please verify that Matlab 2017a or later and the following toolboxes are installed!
 <br /> -- Statistical and machine learning toolbox
 <br /> -- Signal processing toolbox

## Usage of optcorr_batch():

  Run matlab, in command window, type in: (choose one of the following)
  
<br /> **optcorr_batch();**  
call this function without any parameters will popup a FileOpenDialog to select the input file in elmaven output format (.csv) and save the outcome as (_ optcorr.xlsx) in the same folder

<br /> **optcorr_batch(name1, value1, name2, value2, ...);**    
this option allows you to use non-default settings by name-value pairs. For example:

<br />  **optcorr_batch('fname', 'C:\data\example1.csv');**   
if fname is specified, FileOpenDilaglog won't popup.

<br />  **optcorr_batch('resolution', 480000, 'ppm', 3);**    
set instrument resolution and ppm.  default resolution is 140000, default ppm is 5. (make sure the same ppm value is used in elmaven).  

<br />  **optcorr_batch('resolution', 70000, 'purity', 0.98);**  
set instrument resolution and tracer purity.  default purity is 0.99.  No need to specify tracer type, which will be recognized automatically.

<br />  **optcorr_batch('solver', 'optcorr');**  
3 different solvers are available for now, based on different assumptions.   
<br />**'optcorr'** is the default solver, assumes gaussian-shaped peak superposition for unresolved peaks.
<br />**'isocorr'** is the earlier algorithm that does not correct for non-tracer natural abundances and is resolution independant. This is an approximation at infinitly high  resolution.
<br />**'Accucor'** is based on the assumption of simple summation of unresolved peaks (resolution dependant) as the actural observed signal. The results should be similar to Accucor https://github.com/XiaoyangSu/AccuCor and isocorv2 https://github.com/MetaSys-LISBP/IsoCor/


**or write a simple reusable script like this:**
<br /> fname = 'example1.csv';
<br /> res = 480000;
<br /> ppm = 3;
<br /> purity = 0.99;
<br /> solver = 'optcorr';
<br /> optcorr_batch('fname', fname, 'resolution', res, 'ppm', ppm, 'purity', purity, 'solver', solver);
## Input:
Input file is the MAVEN/elMaven peak detection output in .csv format, with isotope detection option enabled.  The isotope for correction is automatically determined by the contents of the isotopeLabel column in your input data. Note that, only single isotope labeled is supported.  see "example1.csv".
## Output: 
The output is an excel file (.xlsx) of the same filename followed by _\_optcorr,_ containing 4 spread sheets:
1. Original: A copy of the original input
2. Corrected: absolute intensities after correction
3. Normalized:  Normalized intensities after correction
4. PoolAfterDF:  sum of all isotopes forms after correction 
5. logs: parameter settings and runtime.

