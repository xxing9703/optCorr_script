# optCorr_script
  optCorr is a simulation & optimization based algorithm for natural isotope abudance correction. 
The state-of-the-art resolution-dependent algorithms like Accucorr and IsoCorv2 determine whether a non-tracer natural isotopic component should be corrected or not based on the MS resolution provided and define the proper correction matrix accordingly to correct only for those "unresolved" ones.  It is assumed that, the observed pre-correction signal is the simple summation of the unresolved peak abundences. However, this is not true since peak top (instead of peak area) is often used as the signal readout when the profile spectrum is turned into centroid.  In fact, the unresolved contaminant peaks only partially contributes to the observed signal, but are 100% subtracted in Accurcorr and IsoCorv2, which cause over correction. optCorr solved this problem from the root by using a built-in MS-iso simulator that simulates the MS profile spectrum of an arbitrary labeling% at a given resolution, from which a centroid spectrum and isotopic distribution signals are extracted and compared to the experimental. the right amount of correction is automatically taken care of by simulation. True labeling% is obtained by non-linear optimiation that minimized the least squre error between the simulated and observed signals. 

<br />  This is a light script verion of optCorr that contains the core "optcorr" function and a wrapper for batch processing elmaven output files: 

<br />   function [distout,err,ysim]=optcorr(mz,atoms,type,abd,imp,fwhm,dmz,distin,option)
<br />   function optcorr_batch(varargin)  

## Usage of optcorr_batch():

  Run matlab, in command window, type 

<br /> **optcorr_batch();**   this will use all default settings and use FileOpenDiaglog to load in elmaven format (.csv) and output the corrected (.xlsx)
<br /> **optcorr_batch(name1,value1,name2,value2,...);**    this option allows you to use non-default settings by name-value pairs. examples:
<br />  **optcorr_batch('fname','C:\data\example1.csv');**   if fname is specified, FileOpenDilaglog won't popup.
<br />  **optcorr_batch('resolution',480000,'ppm',3);**    default resolution is 140000, default ppm is 5.  (use the same ppm setting in elmaven peak picking).  
<br />  **optcorr_batch('resolution',70000,'purity',0.98);**   default purity is 0.99.  no need to specify tracer type, which will be recognized automatically.
<br />  **optcorr_batch('solver','isocorr');**   isocorr is the earlier algorithm that does not correct for non-tracer natural abundances.

 
      

