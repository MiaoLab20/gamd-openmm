# OpenMM - Gaussian Accelerated Molecular Dynamics (GaMD) module

Gaussian Accelerated Molecular Dynamics (GaMD) is a biomolecular enhanced sampling method that works by adding a harmonic boost potential to smoothen the system potential energy surface. By constructing a boost potential that follows Gaussian distribution, accurate reweighting of the GaMD simulations is achieved using cumulant expansion to the second order. GaMD has been demonstrated on three biomolecular model systems: alanine dipeptide, chignolin folding and ligand binding to the T4-lysozyme. Without the need to set predefined reaction coordinates, GaMD enables unconstrained enhanced sampling of these biomolecules. Furthermore, the free energy profiles obtained from reweighting of the GaMD simulations allow us to identify distinct low energy states of the biomolecules and characterize the protein folding and ligand binding pathways quantitatively.


## Installation
This module has been tested with the most recent OpenMM installed using the instructions instructions in the [OpenMM User Guide - Section 2.2 Installing OpenMM](http://docs.openmm.org/latest/userguide/application.html#installing-openmm).

## Status

This code base is currently still in development, has bugs, and is not ready for production use.  You may get incorrect results, while we are still in the development process.  Use at your own risk!
