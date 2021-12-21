# OpenMM - Gaussian Accelerated Molecular Dynamics (GaMD) module

Gaussian Accelerated Molecular Dynamics (GaMD) is a biomolecular enhanced sampling method that works by adding a harmonic boost potential to smoothen the system potential energy surface. By constructing a boost potential that follows Gaussian distribution, accurate reweighting of the GaMD simulations is achieved using cumulant expansion to the second order. GaMD has been demonstrated on three biomolecular model systems: alanine dipeptide, chignolin folding and ligand binding to the T4-lysozyme. Without the need to set predefined reaction coordinates, GaMD enables unconstrained enhanced sampling of these biomolecules. Furthermore, the free energy profiles obtained from reweighting of the GaMD simulations allow us to identify distinct low energy states of the biomolecules and characterize the protein folding and ligand binding pathways quantitatively.


## Installation
1.  You will need to start by installing [Anaconda Python 3.x](https://www.anaconda.com/products/individual#Downloads).
2.  Next, install OpenMM using the instructions found in the [OpenMM User Guide - Section 2.2 Installing OpenMM](http://docs.openmm.org/latest/userguide/application.html#installing-openmm).
3.  You'll need the AmberTools for doing the post MD analysis.  You can do this by executing the following command: 
    ```
    conda install -c conda-forge ambertools=20
    ```  
4.  You'll need the PyReweighting scripts, which can be cloned from 
the [PyReweighting Git Repository](https://github.com/MiaoLab20/PyReweighting).  (NOTE:  If you are 
doing development on the GaMD module itself and want to use the test scripts, the PyReweighting project
 directory should be added to your path, so that the scripts can find it.)
5.  Clone and Install this package: 
    ```
    git clone https://github.com/MiaoLab20/GaMD-OpenMM.git
    cd GaMD-OpenMM
    setup.py install
    ```

## Testing (Optional)
You may also optionally run tests: 
    ```
    setup.py test
    ```

## Run
You can try a test run of GaMD in OpenMM. From within the GaMD-OpenMM/ 
directory:
```
cd gamd
gamdRunner.py xml tests/data/dip_amber.xml
```  

This will made a directory named output/ in the current directory where one
can find all of the GaMD output logs, trajectories, etc. This is a very short
run for demonstration and testing purposes only.

### Important Options and Hints

* The gamdRunner.py program can be run with the '-h' argument to see all
available options. Please see *link to RTD here* for a
detailed description of programs and options.

## Status

     The dihedral and total boost algorithms are working and have been validated.

## Authors and Contributors

The following people have contributed directly to the coding and validation
efforts of GaMD-OpenMM (listed an alphabetical order of last name). 
Thanks also to everyone who has helped or will help improve this project by 
providing feedback, bug reports, or other comments.

* Matthew Copeland
* Hung Do
* Keya Joshi
* Yinglong Miao
* Lane Votapka
* Jinan Wang

### Citing GaMD-OpenMM

If you use GaMD-OpenMM, please cite the following paper:

* CITATION HERE
