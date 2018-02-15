"""
gamd.py: Implements the GaMD integration method.

Portions copyright (c) 2018 University of Kansas
Authors: Matthew Copeland, Yinglong Miao
Contributors:

"""

from __future__ import absolute_import

__author__ = "Matthew Copeland"
__version__ = "1.0"

from simtk.openmm import CustomIntegrator
from simtk.unit import *
from abc import ABCMeta
from abc import abstractmethod


class GamdIntegratorBase(CustomIntegrator):
    __metaclass__ = ABCMeta

    """
        GamdIntegrator implements the GaMD integration algorithm, mode 2 ("boosting the total potential energy only")

        Based on the following papers:
        * J. Chem Theory Comput. 2017, 13, 9-19  - DOI: 10.1021/acs.jctc.6b00931
        * J. Chem Theory Comput. 2015, 11, 3584-3595 - DOI: 10.1021/acs.jctc.5b00436

    """

    # def __init__(self,dt,alpha,E):
    def __init__(self, dt=2.0 * femtoseconds, ntcmd=1000000, nteb=1000000, ntave=10000,
                 sigma0=6.0 * kilocalories_per_mole):
        CustomIntegrator.__init__(self, dt)

        if ntcmd < ntave or ntcmd % ntave != 0:
            raise Exception

        self.dt = dt
        self.ntcmd = ntcmd
        self.nteb = nteb
        self.ntave = ntave
        self.sigma0 = sigma0
        self.step_to_begin_adding_boost_potential = ntcmd
        self.step_to_stop_updating_v_values = ntcmd + nteb

        self.add_common_variables()
        self.addGlobalVariables()
        self.addUpdateContextState()

        # Compute Instructions
        #
        self.addComputeGlobal("currentLocation", "0")
        self.addComputeGlobal("count", "count + 1")
        self.addComputeGlobal("wcount", "wcount + 1")
        self.addComputeGlobal("currentEnergy", "energy")

        # Stage 1
        self.debugStep(1.0)  # 1.0
        self.beginIfBlock("count <= " + str(ntcmd))
        # Conventional / aMD Run
        self.addComputePerDof("v", "v+dt*fprime/m; fprime=f*((1.0-modify) + modify*(alpha/(alpha+E-energy))^2); modify=step(E-energy)")
        self.addComputePerDof("oldx", "x")
        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-oldx)/dt")

        self.endBlock()

        # Stage 1 / 2:

        # * During Stage 1 and 2, this handles calculating the window Potential States
        # * On the ntave window, it does the following:
        #   * Copy the
        #
        #
        self.debugStep(.5) # 1.5
        self.beginIfBlock("count <= " + str(self.step_to_stop_updating_v_values))

        self.addEnergyValueCalculations()
        self.beginIfBlock("wcount >= " + str(ntave))

        self.update_potential_state_values_with_window_potential_state_values()
        self.add_calculate_E_k0()

        self.endBlock()  # wcount >= ntave
        self.endBlock()  # count <= step_to_stop_updating_v_values




        # Stage 2 / 3
        self.debugStep(.5)  # 2.0
        self.beginIfBlock("count >= " + str(self.step_to_begin_adding_boost_potential))

        self.addComputePerDof("newx", "0.0")
        self.addComputePerDof("newv", "0.0")
        self.addComputePerDof("oldx", "0.0")
        self.addComputePerDof("oldv", "0.0")

        # Save off Debugging Values
        self.addComputePerDof("forceOfF", "f")
        self.addComputePerDof("mass", "m")

        #
        # Do GaMD Steps to calculate the boostPotential and the scaling Force
        # and then update the position and velocity
        #
        self.debugStep(.01)  # 2.01
        self.addComputeGlobal("boostPotential", "0.5 * k0 * (E-energy)^2/(Vmax-Vmin)")
        self.debugStep(.01)  # 2.02
        self.addComputeGlobal("boostForce", "1.0-((k0* (E - energy))/(Vmax - Vmin)) ")
        self.debugStep(.01)  # 2.03
        self.addComputePerDof("scalingForce", "f*boostForce")
        self.debugStep(.01)  # 2.04
        self.addComputePerDof("scalingVelocity", "v + scalingForce*dt/m")
        self.debugStep(.01)  # 2.05
        self.addComputePerDof("oldv", "v")
        self.debugStep(.01)  # 2.06
        self.addComputePerDof("oldx", "x")
        self.debugStep(.01)  # 2.07
        self.addComputePerDof("newx", "x + dt * scalingVelocity")
        self.debugStep(.01)  # 2.08
        self.addComputePerDof("newv", "(newx - oldx)/dt")
        self.debugStep(.01)  # 2.09
        self.addComputePerDof("x","newx")
        self.debugStep(.01)  # 2.10
        self.addComputePerDof("v","newv")
        self.debugStep(.01)  # 2.11
        # Save off Debugging Values

        self.endBlock()
        self.addConstrainPositions()
        self.addComputePerDof("XafterContrainPosition", "x")
        self.debugStep(.01)   # 2.12


    def debugStep(self, number):
        self.addComputeGlobal("currentLocation", "currentLocation + " + str(number))

    def update_potential_state_values_with_window_potential_state_values(self):
        # Update window variables
        self.addComputeGlobal("Vmax", "wVmax")
        self.addComputeGlobal("Vmin", "wVmin")
        self.addComputeGlobal("Vavg", "wVavg")
        self.addComputeGlobal("sigmaV", "sqrt(wVariance)")

        # reset variables
        self.addComputeGlobal("wcount", "0")
        self.addComputeGlobal("M2", "0")
        self.addComputeGlobal("wVmax", "-1E99")
        self.addComputeGlobal("wVmin", "1E99")
        self.addComputeGlobal("wVavg", "0")
        self.addComputeGlobal("oldVavg", "0")
        self.addComputeGlobal("wVariance", "0")

    def add_common_variables(self):
        # Global Variable Declarations
        self.addGlobalVariable("sigma0", self.sigma0)
        self.addGlobalVariable("alpha", 0)
        self.addGlobalVariable("E", -1E99)
        self.addGlobalVariable("count", 0)
        self.addGlobalVariable("wcount", 0)

        self.addGlobalVariable("Vmax", -1E99)
        self.addGlobalVariable("Vmin", 1E99)
        self.addGlobalVariable("Vavg", 0)
        self.addGlobalVariable("oldVavg", 0)
        self.addGlobalVariable("sigmaV", 0)

        self.addGlobalVariable("M2", 0)
        self.addGlobalVariable("wVmax", -1E99)
        self.addGlobalVariable("wVmin", 1E99)
        self.addGlobalVariable("wVavg", 0)
        self.addGlobalVariable("wVariance", 0)

        self.addGlobalVariable("boostPotential", 0)
        self.addGlobalVariable("boostForce", 0)

        self.addPerDofVariable("scalingForce", 0)
        self.addPerDofVariable("scalingVelocity", 0)
        self.addPerDofVariable("oldx", 0)

        # Debugging Values
        self.addGlobalVariable("currentLocation", 0)
        self.addPerDofVariable("oldv", 0)
        self.addPerDofVariable("newv", 0)
        self.addPerDofVariable("newx", 0)
        self.addGlobalVariable("currentEnergy", 0)
        self.addPerDofVariable("forceOfF", 0)
        self.addPerDofVariable("mass", 0)
        self.addPerDofVariable("XafterContrainPosition",0)

    @abstractmethod
    def addGlobalVariables(self):
        raise NotImplementedError("must implement addGlobalVariables")

    @abstractmethod
    def addEnergyValueCalculations(self):
        raise NotImplementedError("must implement addEnergyValueCalculations")

    @abstractmethod
    def add_calculate_E_k0(self):
        raise NotImplementedError("must implement calculate_E_k0")

    @abstractmethod
    def printVariables(self):
        raise NotImplementedError("must implement printVariables")


class GamdIntegratorBoostTotalPotential(GamdIntegratorBase):

    def __init__(self, dt, ntcmd, nteb, ntave, sigma0):
        GamdIntegratorBase.__init__(self, dt, ntcmd, nteb, ntave, sigma0)

    def addGlobalVariables(self):
        self.addGlobalVariable("k0", 0)
        self.addGlobalVariable("k0prime", 0)
        self.addGlobalVariable("k0doubleprime", 0)

    def addEnergyValueCalculations(self):
        self.addComputeGlobal("wVmax", "max(energy, wVmax)")
        self.addComputeGlobal("wVmin", "min(energy, wVmin)")

        self.addComputeGlobal("oldVavg", " wVavg")
        self.addComputeGlobal("wVavg", "wVavg + (energy-wVavg)/wcount")
        self.addComputeGlobal("M2", "M2 + (energy - oldVavg)*(energy - wVavg)")
        self.addComputeGlobal("wVariance", "select(wcount-1,M2/(wcount - 1),0)")

    @abstractmethod
    def add_calculate_E_k0(self):
        raise NotImplementedError("must implement calculate_E_k0")

    def printVariables(self):
        wvmax = str(self.getGlobalVariableByName("wVmax"))
        wvmin = str(self.getGlobalVariableByName("wVmin"))
        wvavg = str(self.getGlobalVariableByName("wVavg"))
        wvariance = str(self.getGlobalVariableByName("wVariance"))
        sigmav = str(self.getGlobalVariableByName("sigmaV"))
        vmax = str(self.getGlobalVariableByName("Vmax"))
        vmin = str(self.getGlobalVariableByName("Vmin"))
        vavg = str(self.getGlobalVariableByName("Vavg"))

        m2 = str(self.getGlobalVariableByName("M2"))
        k0 = str(self.getGlobalVariableByName("k0"))
        k0prime = str(self.getGlobalVariableByName("k0prime"))
        E = str(self.getGlobalVariableByName("E"))
        count = str(self.getGlobalVariableByName("count"))
        wcount = str(self.getGlobalVariableByName("wcount"))
        boostPotential = str(self.getGlobalVariableByName("boostPotential"))
        boostForce = str(self.getGlobalVariableByName("boostForce"))

        print("\nTime Window " + str(count) + ":")
        print("-------------------------------")

        print("\nWindow Variables:  ")
        print("count, wcount, wVmax,         wVmin,           wVavg,         wVariance,         M2")
        print(count + ", " + wcount + ", " + wvmax + ", " + wvmin + ", " + wvavg + ", " + wvariance + ", " + m2)

        print("\nCalculated Window Variables:")
        print("Vmax,              Vmin,         Vavg,         sigmaV,         k0,      k0prime,         E")
        print(vmax + ", " + vmin + ", " + vavg + ", " + sigmav + ", " + k0 + ", " + k0prime + ", " + E)

        print("\nCalculation Variables:")
        print("current potential, boostPotential, boostForce")
        print(str(self.getGlobalVariableByName("currentEnergy")) + ",    " + boostPotential + ",  " + boostForce)

    def printPositions(self):
        count = str(self.getGlobalVariableByName("count"))

        print("\nCurrent Debug Location:  " + str(self.getGlobalVariableByName("currentLocation")))
        print("Position Variables + " + count + ":")
        self.printDofVariable("oldx")
        self.printDofVariable("newx")
        self.printDofVariable("XafterContrainPosition")

        #self.printDofVariable("forceOfF")
        #self.printDofVariable("mass")
        #self.printDofVariable("oldv")
        #self.printDofVariable("newv")
        #self.printDofVariable("scalingForce")
        #self.printDofVariable("scalingVelocity")

    def printDofVariable(self, name):
        dofValue = self.getPerDofVariableByName(name)
        print("\n" + name + ":\n")
        for i in range(len(dofValue)):
            print(dofValue[i])


class GamdIntegratorBoostTotalPotentialLowerBound(GamdIntegratorBoostTotalPotential):

    def __init__(self, dt, ntcmd, nteb, ntave, sigma0):
        GamdIntegratorBoostTotalPotential.__init__(self, dt, ntcmd, nteb, ntave, sigma0)

    def add_calculate_E_k0(self):
        self.addComputeGlobal("E", "Vmax")
        self.addComputeGlobal("k0prime", "(sigma0/sigmaV) * (Vmax - Vmin)/(Vmax - Vavg)")
        self.addComputeGlobal("k0", "min(1.0, k0prime); ")


class GamdIntegratorBoostTotalPotentialUpperBound(GamdIntegratorBoostTotalPotential):

    def __init__(self, dt, ntcmd, nteb, ntave, sigma0):
        GamdIntegratorBoostTotalPotential.__init__(self, dt, ntcmd, nteb, ntave, sigma0)

    def add_calculate_E_k0(self):
        self.addComputeGlobal("k0", "1.0")  # Our Else Case for the If Block
        self.addComputeGlobal("k0doubleprime", "(1-(sigma0/sigmaV) * (Vmax - Vmin)/(Vavg - Vmin)")
        self.beginIfBlock("0.0 < k0doubleprime <= 1.0")
        self.addComputeGlobal("k0", "k0doubleprime")
        self.addComputeGlobal("E", "Vmin + (Vmax - Vmin)/k0")
        self.endBlock()
