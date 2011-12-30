== Layout of src code ==
 * /lib/src   - subdirectory for each project (hold .c, .h, and individual makefiles)
 * /lib/build - hold build scripts for generating dll's and so's 
 * /lib       - hold xx.py wrapper code and the dll/so's

== Build procedure MSVS ==
 * You need the following environment variables:
    - MSVSHOME : path to MS Visual Studio 
    - TDLDIR : /path/to/tdl
 * cd into lib/build
 * exectute the build script: build_msvs.bat 

== Build procedure *nix ==
 * cd into lib/build 
 * run make  

== Modules and dependencies ==
  - get rid of all extern dependencies (put gsl stuff in numfcns.c)  
 * GSL
 * XRR: Ifeffit.py (http://cars9.uchicago.edu/ifeffit)


