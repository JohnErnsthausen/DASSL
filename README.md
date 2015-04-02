# DASSL 
### ODE and DAE Numerical Integration Software

Linda Petzold wrote DASSL.

CWI testset drivers and documentation are available from their website.

[CWI](http://www.dm.uniba.it/~testset)

I rewrote a few drivers for ODE and DAE examples to familiarize myself with DASSL. A purpose
of this repository is to share the drivers. I also share a rakefile to demonstrate building
a DASSL library and building a DASSL driver and the availability of a fortran 77 software development
environment.


The gfortran compiler in the package

     gcc version 4.6.3 (Ubuntu/Linaro 4.6.3-1ubuntu5)

will compile all the codes in this repository. I was surprised too! I thought every numerical
code since 2000 had to be written in C++ or Python or MatLab. 

I assume your system has Ruby and Rake installed. These references will walk you through installing ruby and rake.

* [Ruby installation instructions](https://www.ruby-lang.org)
* [Rake installation instructions](http://rake.rubyforge.org)


In your dassl directory, change to the tests directory.

This command will list all the rake commands available to you.

     $ rake -T

Then compile the dassl library via

     $ rake make

Then compile a driver

     $ rake compile[chemakzo.f]

Run the examples as usual

     $ ./chemakzo.exe

To remove .dat files and .o files and .exe files

     $ rake clean

To remove everything a clean command will remove and restore the repository to its original state, run

     $ rake clobber

The run and datamine tasks will run the example on a mesh of tolerances and datamine the results.

Enjoy! Email me if you like the software or have a comment.

John Ernsthausen
