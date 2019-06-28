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
            raise ValueError("ntcmd must be greater than and a multiple of ntave.")

        self.dt = dt
        self.ntcmd = ntcmd
        self.nteb = nteb
        self.ntave = ntave
        self.sigma0 = sigma0
        self.step_to_begin_adding_boost_potential = ntcmd
        self.step_to_stop_updating_v_values = ntcmd + nteb

        #
        # Variable Definitions
        #
        self.add_common_variables()
        self.add_global_variables()
        self.addUpdateContextState()

        #
        # Common Compute Instructions
        #
        self.addComputeGlobal("count", "count + 1")
        self.addComputeGlobal("wcount", "wcount + 1")
        self.addComputeGlobal("currentEnergy", "energy")
        self.addComputeGlobal("Vmax", "max(Vmax,energy)")
        self.addComputeGlobal("Vmin", "min(Vmin,energy)")

        self.add_stage_1_conventional_md_instructions()
        self.add_common_stage_1_and_2_instructions()
        self.add_stage_2_and_3_instructions()

        self.addComputePerDof("coordinates","x")

    def add_common_variables(self):
        # Global Variable Declarations
        self.addGlobalVariable("sigma0", self.sigma0)
        self.addGlobalVariable("alpha", 0)
        self.addGlobalVariable("E", -1E99)
        self.addGlobalVariable("count", 0)
        self.addGlobalVariable("wcount", 0)
        self.addGlobalVariable("k0", 0)

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
        self.addGlobalVariable("totalForceScalingFactor", 0)

        self.addPerDofVariable("scalingForce", 0)
        self.addPerDofVariable("scalingVelocity", 0)
        self.addPerDofVariable("oldx", 0)
        self.addPerDofVariable("coordinates",0)

        # Debugging Values
        self.addPerDofVariable("oldv", 0)
        self.addPerDofVariable("newv", 0)
        self.addPerDofVariable("newx", 0)
        self.addGlobalVariable("currentEnergy", 0)

    def add_stage_1_conventional_md_instructions(self, group=''):
        # Stage 1
        self.beginIfBlock("count <= " + str(self.step_to_begin_adding_boost_potential))

        # Conventional / aMD Run
        self.addComputePerDof("fg", "f"+str(group))
        self.addComputePerDof("v", "v+dt*fprime/m; fprime=fg*((1.0-modify) + modify*(alpha/(alpha+E-energy))^2); modify=step(E-energy)" % group)
        self.addComputePerDof("oldx", "x")
        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-oldx)/dt")

        self.endBlock()

    def add_common_stage_1_and_2_instructions(self):
        # Stage 1 / 2:
        # * During Stage 1 and 2, this handles calculating the window Potential States
        # * On the ntave window, it does the following:
        #   * Copy the
        #
        #
        self.beginIfBlock("count <= " + str(self.step_to_stop_updating_v_values))
        self.add_energy_value_calculations()
        self.beginIfBlock("wcount >= " + str(self.ntave))
        self.update_potential_state_values_with_window_potential_state_values()
        #
        #
        # If we start modifying E prior to the aMD calculations completion, then
        # it will blow up.
        #
        self.beginIfBlock("count >=" + str(self.step_to_begin_adding_boost_potential))
        self.add_calculate_E_k0()
        self.endBlock()
        self.endBlock()  # wcount >= ntave
        self.endBlock()  # count <= step_to_stop_updating_v_values

    def update_potential_state_values_with_window_potential_state_values(self):
        # Update window variables
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

    def add_stage_2_and_3_instructions(self):
        #
        # Stage 2 / 3
        #
        self.beginIfBlock("count >= " + str(self.step_to_begin_adding_boost_potential))
        self.addComputePerDof("newx", "0.0")
        self.addComputePerDof("newv", "0.0")
        self.addComputePerDof("oldx", "0.0")
        self.addComputePerDof("oldv", "0.0")
        #
        # Do GaMD Steps to calculate the boostPotential and the scaling Force
        # and then update the position and velocity
        #
        self.addComputeGlobal("boostPotential", "0.5 * k0 * (E-energy)^2/(Vmax-Vmin)")
        self.addComputeGlobal("totalForceScalingFactor", "1.0-((k0 * (E - energy))/(Vmax - Vmin)) ")
        self.addComputePerDof("scalingForce", "f*totalForceScalingFactor")
        self.addComputePerDof("scalingVelocity", "v + scalingForce*dt/m")
        self.addComputePerDof("oldv", "v")
        self.addComputePerDof("oldx", "x")
        self.addComputePerDof("newx", "x + dt * scalingVelocity")
        self.addComputePerDof("x", "newx")
        self.addConstrainPositions()
        self.addComputePerDof("newv", "(newx - oldx)/dt")
        self.addComputePerDof("v", "newv")
        self.addConstrainVelocities()
        self.endBlock()  # count >= step_to_begin_adding_boost_potential

    def get_boost_potential(self):
        return self.getGlobalVariableByName("boostPotential")

    def get_current_potential_energy(self):
        return self.getGlobalVariableByName("currentEnergy")

    def get_total_force_scaling_factor(self):
        return self.getGlobalVariableByName("totalForceScalingFactor")

    def get_coordinates(self):
        return self.getPerDofVariableByName("coordinates")

    @abstractmethod
    def add_global_variables(self):
        raise NotImplementedError("must implement addGlobalVariables")

    @abstractmethod
    def add_energy_value_calculations(self):
        raise NotImplementedError("must implement addEnergyValueCalculations")

    @abstractmethod
    def add_calculate_E_k0(self):
        raise NotImplementedError("must implement calculate_E_k0")


class GamdTotalBoostPotentialIntegrator(GamdIntegratorBase):

    def __init__(self, dt, ntcmd, nteb, ntave, sigma0):
        GamdIntegratorBase.__init__(self, dt, ntcmd, nteb, ntave, sigma0)

    def add_global_variables(self):

        self.addGlobalVariable("k0prime", 0)
        self.addGlobalVariable("k0doubleprime", 0)

    def add_energy_value_calculations(self):
        self.addComputeGlobal("wVmax", "max(energy, wVmax)")
        self.addComputeGlobal("wVmin", "min(energy, wVmin)")

        self.addComputeGlobal("oldVavg", " wVavg")
        self.addComputeGlobal("wVavg", "wVavg + (energy-wVavg)/wcount")
        self.addComputeGlobal("M2", "M2 + (energy - oldVavg)*(energy - wVavg)")
        self.addComputeGlobal("wVariance", "select(wcount-1,M2/(wcount - 1),0)")

    @abstractmethod
    def add_calculate_E_k0(self):
        raise NotImplementedError("must implement calculate_E_k0")


class GamdTotalBoostPotentialIntegratorLowerBound(GamdTotalBoostPotentialIntegrator):

    def __init__(self, dt, ntcmd, nteb, ntave, sigma0):
        GamdTotalBoostPotentialIntegrator.__init__(self, dt, ntcmd, nteb, ntave, sigma0)

    def add_calculate_E_k0(self):
        self.addComputeGlobal("E", "Vmax")
        self.addComputeGlobal("k0prime", "(sigma0/sigmaV) * (Vmax - Vmin)/(Vmax - Vavg)")
        self.addComputeGlobal("k0", "min(1.0, k0prime); ")


class GamdTotalBoostPotentialIntegratorUpperBound(GamdTotalBoostPotentialIntegrator):

    def __init__(self, dt, ntcmd, nteb, ntave, sigma0):
        GamdTotalBoostPotentialIntegrator.__init__(self, dt, ntcmd, nteb, ntave, sigma0)

    def add_calculate_E_k0(self):
        self.addComputeGlobal("k0", "1.0")  # Our Else Case for the If Block
        self.addComputeGlobal("k0doubleprime", "1-(sigma0/sigmaV) * (Vmax - Vmin)/(Vavg - Vmin)")
        self.beginIfBlock("0.0 < k0doubleprime <= 1.0")
        self.addComputeGlobal("k0", "k0doubleprime")
        self.addComputeGlobal("E", "Vmin + (Vmax - Vmin)/k0")
        self.endBlock()
        
