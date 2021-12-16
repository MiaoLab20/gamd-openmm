.. gamd documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GaMD-OpenMM: Gaussian Accelerated Molecular Dynamics using OpenMM
=================================================================

Gaussian Accelerated Molecular Dynamics (GaMD) is a biomolecular enhanced 
sampling method that works by adding a harmonic boost potential to smoothen 
the system potential energy surface. By constructing a boost potential that 
follows Gaussian distribution, accurate reweighting of the GaMD simulations 
is achieved using cumulant expansion to the second order. GaMD has been 
demonstrated on three biomolecular model systems: alanine dipeptide, 
chignolin folding and ligand binding to the T4-lysozyme. Without the need 
to set predefined reaction coordinates, GaMD enables unconstrained enhanced 
sampling of these biomolecules. Furthermore, the free energy profiles obtained 
from reweighting of the GaMD simulations allow us to identify distinct low 
energy states of the biomolecules and characterize the protein folding and 
ligand binding pathways quantitatively.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   running
   input_files
   tutorial
   api

Cite GaMD-OpenMM
================

If you wish to cite GaMD-OpenMM, please cite the following paper:

* GAMD OPENMM PAPER

You may also optionally cite one or more of the following papers:

* MORE GAMD PAPERS

Getting Involved
================

Please report **bugs** or **enhancement requests** through the `Issue 
Tracker`_. 

.. _Issue Tracker: https://github.com/MiaoLab20/GaMD-OpenMM/issues

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
