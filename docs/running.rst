Running GaMD-OpenMM
===================

You can try a test run of GaMD in OpenMM. From within the GaMD-OpenMM/ 
directory::

  cd gamd
  python gamdRunner.py xml tests/data/dip_amber.xml

This will made a directory named output/ in the current directory where one
can find all of the GaMD output logs, trajectories, etc. This is a very short
run for demonstration and testing purposes only. You will need to make and 
modify your own XML input file for your system of interest. Please see
the :doc:`Input files<input_files>` documentation for more information and
guidance.

The gamdRunner.py program can be run with the '-h' argument to see all
available options.
