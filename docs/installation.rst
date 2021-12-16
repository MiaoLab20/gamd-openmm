Installation
============

At this time, GaMD-OpenMM has only been tested on Linux systems. Therefore, all
installation instructions are for Linux only.

Install Conda
-------------

It is recommended, though not mandatory, that you install Conda to use with 
GaMD-OpenMM. Without Conda, all dependencies will need to be installed by hand.

If you do not already have Conda, it can be easily installed by completing the
following steps:

Download Conda with Python version 3.8 from the website. 
https://www.anaconda.com/products/individual#Downloads. Alternatively, if you
wish you use the command line, you can execute the following commands and 
fill out the prompts::

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  sh Miniconda3-latest-Linux-x86_64.sh

Make sure Conda is installed by running:

``which conda``

You will want to use Python 3.8, so you can see which version you are with
the command:

``python -V``

If it says any other version besides Python 3.8, then enter:

``conda install python=3.8``

If you want you can create a conda environment, 

``conda create --name GAMD python=3.8``

but you can also just install all packages straight to the base environment
if you wish to. If using an environment, whenever you're installing or running 
anything involving OpenMM or GaMD-OpenMM, make sure that you have activated your 
environment by running ``conda activate GAMD``.

Install OpenMM with Conda
------------------------------------
This section describes the fasted and easiest way to get GaMD-OpenMM working.

You must install OpenMM either from conda or from source. If you wish to 
install from source, see the documentation at 
[OpenMM User Guide - Section 2.2 Installing OpenMM](http://docs.openmm.org/latest/userguide/application.html#installing-openmm).

With Conda working, You may create and activate any environment you wish, 
or use the base environment. Install OpenMM:

``conda install -c conda-forge openmm``

Installation of GaMD-OpenMM itself begins with cloning and installing the 
GaMD-OpenMM python API::

  git clone https://github.com/MiaoLab20/GaMD-OpenMM.git
  cd GaMD-OpenMM
  python setup.py install

Once OpenMM and the GaMD-OpenMM are installed, it is recommended that 
you run tests of GaMD-OpenMM. From within the "GaMD-OpenMM/" directory, run:

``python setup.py test``

You should now be able to use GaMD-OpenMM.