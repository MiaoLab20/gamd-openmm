# OpenMM - Gaussian Accelerated Molecular Dynamics (GaMD) module

Gaussian Accelerated Molecular Dynamics (GaMD) is a biomolecular enhanced sampling method that works by adding a harmonic boost potential to smoothen the system potential energy surface. By constructing a boost potential that follows Gaussian distribution, accurate reweighting of the GaMD simulations is achieved using cumulant expansion to the second order. GaMD has been demonstrated on three biomolecular model systems: alanine dipeptide, chignolin folding and ligand binding to the T4-lysozyme. Without the need to set predefined reaction coordinates, GaMD enables unconstrained enhanced sampling of these biomolecules. Furthermore, the free energy profiles obtained from reweighting of the GaMD simulations allow us to identify distinct low energy states of the biomolecules and characterize the protein folding and ligand binding pathways quantitatively.


## Installation
1.  You will need to start by installing [Anaconda Python 3.x](https://www.anaconda.com/products/individual#Downloads).
2.  Next, install OpenMM using the instructions found in the [OpenMM User Guide - Section 2.2 Installing OpenMM](http://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm).
3.  You'll need the AmberTools for doing the post MD analysis.  You can get these tools by executing the following command: 
    ```
    conda install -c conda-forge ambertools
    ```  
4.  You'll need the PyReweighting scripts, which can be cloned from 
the [PyReweighting Git Repository](https://github.com/MiaoLab20/PyReweighting).  (NOTE:  If you are 
doing development on the GaMD module itself and want to use the test scripts, the PyReweighting project
 directory should be added to your path, so that the scripts can find it.)
5.  Clone and Install this package: 
    ```
    git clone https://github.com/MiaoLab20/gamd-openmm.git
    cd gamd-openmm
    setup.py install
    ```
6.  The command gamdRunner can either be copied into your user bin directory or you can updated your
PATH variable to include the location fo the gamd-openmm directory, if you would like to use the
gamdRunner for running your simulations.


## Testing (Optional)
You may also optionally run tests: 
    ```
    setup.py test
    ```

## Run

You can run gamd by providing your own configuration file to the gamdRunner
program like the example here.

```
gamdRunner xml configuration-file.xml
```  

We've created the repository [gamd-openmm-examples](https://github.com/MiaoLab20/gamd-openmm-examples) which
contains examples (data files and configuration files) and instructions you can use to validate
your gamd installation. This project can also help you learn how to use the available command line
options to the gamdRunner and about some of the options for the configuration file.

*NOTE:*  The gamdRunner and the gamd-openmm code currently only supports running the conventional md, equilibration, and production stages as a part of a single execution.

### Important Options and Hints

* The gamdRunner program can be run with the '-h' argument to see all
available options. Please see *link to RTD here* for a
detailed description of programs and options.

## Status

We have implemented the upper and lower bound versions of the following types of
gamd boosts:

* dihedral
* total
* dual total/dihedral
* non-bonded
* dual non-bonded/dihedral


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

* <ins>Copeland, M., Do, HN,</ins> Votapka, L., Joshi, K., Wang, J., Amaro, R., and Miao, Y.* (2022) Gaussian accelerated molecular dynamics in OpenMM. The Journal of Physical Chemistry B. In Review. 
