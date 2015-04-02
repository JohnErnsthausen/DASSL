DASSL DAE solver software and Drivers
=====================================

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
