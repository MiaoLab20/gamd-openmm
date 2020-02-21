"""
gamd.py: Implements the GaMD integration method.

Portions copyright (c) 2018 University of Kansas
Authors: Matthew Copeland, Yinglong Miao
Contributors: Lane Votapka

"""

from __future__ import absolute_import

__author__ = "Matthew Copeland"
__version__ = "1.0"

from simtk.openmm import CustomIntegrator
from simtk.unit import *
from abc import ABCMeta
from abc import abstractmethod
import numpy


class GamdIntegratorBase(CustomIntegrator):
    __metaclass__ = ABCMeta

    """
        GamdIntegrator implements the GaMD integration algorithm, all modes

        Based on the following papers:
        * J. Chem Theory Comput. 2017, 13, 9-19  - DOI: 10.1021/acs.jctc.6b00931
        * J. Chem Theory Comput. 2015, 11, 3584-3595 - DOI: 10.1021/acs.jctc.5b00436

    """

    # def __init__(self,dt,alpha,E):
    def __init__(self, dt=2.0 * femtoseconds, ntcmd=1000000, nteb=1000000, ntave=50000,
                 sigma0=6.0 * kilocalories_per_mole):
        """
        Params:
        ntcmd: Number of initial conventional MD steps used to calculate the maximum, minimum, average,
            and standard deviation of the system potential energies. Default is one million timesteps
            for a 2fs timestep.
        nteb: Number of steps to equilibrate the system after adding the boost potential. Default is 
            one million timesteps for a 2fs timestep.
        ntave: Stride to use for averaging for the calculation of Vmax, Vmin, etc. during conventional MD
        """
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

        #
        # Common Compute Instructions
        #
        self.addComputeGlobal("count", "count + 1")
        self.addComputeGlobal("wcount", "wcount + 1")
        
        

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
        
        self.addPerDofVariable("scalingVelocity", 0)
        self.addPerDofVariable("oldx", 0)
        self.addPerDofVariable("coordinates",0)

        # Debugging Values
        self.addPerDofVariable("oldv", 0)
        self.addPerDofVariable("newv", 0)
        self.addPerDofVariable("newx", 0)
        self.addGlobalVariable("currentEnergy", 0)
        self.addGlobalVariable("totalEnergy", 0)
    
    @abstractmethod
    def add_stage_1_conventional_md_instructions(self):
        raise NotImplementedError("add_stage_1_conventional_md_instructions")
    
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
        '''
        self.addComputePerDof("newx", "0.0")
        self.addComputePerDof("newv", "0.0")
        self.addComputePerDof("oldx", "0.0")
        self.addComputePerDof("oldv", "0.0")
        '''

    def get_boost_potential(self):
        return self.getGlobalVariableByName("boostPotential")

    def get_current_potential_energy(self):
        return self.getGlobalVariableByName("totalEnergy")

    def get_coordinates(self):
        return self.getPerDofVariableByName("coordinates")

    def add_global_variables(self):

        self.addGlobalVariable("k0prime", 0)
        self.addGlobalVariable("k0doubleprime", 0)

    def add_energy_value_calculations(self):
        self.addComputeGlobal("wVmax", "max(currentEnergy, wVmax)")
        self.addComputeGlobal("wVmin", "min(currentEnergy, wVmin)")

        self.addComputeGlobal("oldVavg", " wVavg")
        self.addComputeGlobal("wVavg", "wVavg + (currentEnergy-wVavg)/wcount")
        self.addComputeGlobal("M2", "M2 + (currentEnergy - oldVavg)*(currentEnergy - wVavg)")
        self.addComputeGlobal("wVariance", "select(wcount-1,M2/(wcount - 1),0)")
    
    @abstractmethod
    def add_calculate_E_k0(self):
        raise NotImplementedError("must implement calculate_E_k0")


class GamdTotalBoostPotentialLangevinIntegrator(GamdIntegratorBase):

    def __init__(self, dt, ntcmd, nteb, ntave, sigma0, temperature, gamma):
        super(GamdTotalBoostPotentialLangevinIntegrator, self).__init__(dt, ntcmd, nteb, ntave, sigma0)
        kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA
        a = numpy.exp(-gamma*dt) # current velocity component
        b = numpy.sqrt(1 - numpy.exp(- 2 * gamma * dt)) # random velocity component
        #b = numpy.sqrt(kt*(1-gamma*gamma))
        self.kT = kB * temperature
        self.a = a
        self.b = b
        self.gamma = gamma
        self.addGlobalVariable("kT", self.kT) # thermal energy
        self.addGlobalVariable("a", self.a)
        self.addGlobalVariable("b", self.b)
        self.addGlobalVariable("gamma", self.gamma)
        self.addPerDofVariable("sigma", 0) 
        
        
        # currentEnergy is the total energy
        self.addComputeGlobal("currentEnergy", "energy")
        self.addComputeGlobal("totalEnergy", "currentEnergy")
        self.addComputeGlobal("Vmax", "max(Vmax,currentEnergy)")
        self.addComputeGlobal("Vmin", "min(Vmin,currentEnergy)")
        
        self.addUpdateContextState()
        self.add_stage_1_conventional_md_instructions()
        self.add_common_stage_1_and_2_instructions()
        self.add_stage_2_and_3_instructions()
        
        self.addComputePerDof("coordinates","x")
        
    def add_common_variables(self):
        super(GamdTotalBoostPotentialLangevinIntegrator, self).add_common_variables()
        self.addGlobalVariable("totalForceScalingFactor", 0)
        
    def add_stage_1_conventional_md_instructions(self):
        # Stage 1
        self.beginIfBlock("count <= " + str(self.step_to_begin_adding_boost_potential))
        
        # Conventional / aMD Run
        self.addComputePerDof("sigma", "sqrt(kT/m)")
        
        # V step
        self.addComputePerDof("v", "v + 0.5*dt*fprime/m;fprime=f*((1.0-modify) + modify*(alpha/(alpha+E-currentEnergy))^2); modify=step(E-currentEnergy)")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()

        # O step
        self.addComputePerDof("v", "a*v + b*sigma*gaussian")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()
        
        self.addComputeGlobal("currentEnergy", "energy")
        
        # V step
        self.addComputePerDof("v", "v + 0.5*dt*fprime/m;fprime=f*((1.0-modify) + modify*(alpha/(alpha+E-currentEnergy))^2); modify=step(E-currentEnergy)")
        self.addConstrainVelocities()
        
        self.endBlock()
        
    def add_stage_2_and_3_instructions(self):
        super(GamdTotalBoostPotentialLangevinIntegrator, self).add_stage_2_and_3_instructions()
        #
        # Do GaMD Steps to calculate the boostPotential and the scaling Force
        # and then update the position and velocity
        #
        self.addUpdateContextState()
        
        self.addComputeGlobal("boostPotential", "0.5 * k0 * (E-currentEnergy)^2/(Vmax-Vmin)")
        self.addComputeGlobal("totalForceScalingFactor", "1.0-((k0 * (E - currentEnergy))/(Vmax - Vmin)) ")
        
        self.addComputePerDof("sigma", "sqrt(kT/m)")
        
        # V step
        self.addComputePerDof("v", "v + 0.5*(dt*f*totalForceScalingFactor/m)")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()

        # O step
        self.addComputePerDof("v", "a*v + b*sigma*gaussian")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()
        
        self.addComputeGlobal("currentEnergy", "energy")
        self.addComputeGlobal("Vmax", "max(Vmax,currentEnergy)")
        self.addComputeGlobal("Vmin", "min(Vmin,currentEnergy)")
        self.addComputeGlobal("totalForceScalingFactor", "1.0-((k0 * (E - currentEnergy))/(Vmax - Vmin)) ")
        
        # V step
        self.addComputePerDof("v", "v + 0.5*(dt*f*totalForceScalingFactor/m)")
        self.addConstrainVelocities()
        
        
        self.endBlock()  # count >= step_to_begin_adding_boost_potential
        
    def get_total_force_scaling_factor(self):
        return self.getGlobalVariableByName("totalForceScalingFactor")

class GamdGroupBoostPotentialLangevinIntegrator(GamdIntegratorBase):

    def __init__(self, dt, ntcmd, nteb, ntave, sigma0, temperature, gamma, group):
        self.group = group
        super(GamdGroupBoostPotentialLangevinIntegrator, self).__init__(dt, ntcmd, nteb, ntave, sigma0)
        kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA
        a = numpy.exp(-gamma*dt) # current velocity component
        b = numpy.sqrt(1 - numpy.exp(- 2 * gamma * dt)) # random velocity component
        #b = numpy.sqrt(kt*(1-gamma*gamma))
        self.kT = kB * temperature
        self.a = a
        self.b = b
        self.gamma = gamma
        self.addGlobalVariable("kT", self.kT) # thermal energy
        self.addGlobalVariable("a", self.a)
        self.addGlobalVariable("b", self.b)
        self.addGlobalVariable("gamma", self.gamma)
        self.addPerDofVariable("sigma", 0) 
        
        
        # currentEnergy is the total energy
        self.addComputeGlobal("groupEnergy", "energy"+str(group))
        self.addComputeGlobal("currentEnergy", "groupEnergy")
        self.addComputeGlobal("totalEnergy", "energy")
        self.addComputeGlobal("Vmax", "max(Vmax,groupEnergy)")
        self.addComputeGlobal("Vmin", "min(Vmin,groupEnergy)")
        
        self.addUpdateContextState()
        self.add_stage_1_conventional_md_instructions()
        self.add_common_stage_1_and_2_instructions()
        self.add_stage_2_and_3_instructions()
        
        

        self.addComputePerDof("coordinates","x")
        
    def add_common_variables(self):
        super(GamdGroupBoostPotentialLangevinIntegrator, self).add_common_variables()
        self.addGlobalVariable("groupForceScalingFactor", 0)
        self.addGlobalVariable("groupEnergy", 0)
        #self.addGlobalVariable("totalEnergy", 0)
        self.addPerDofVariable("fg", 0)
        
    def add_stage_1_conventional_md_instructions(self):
        # Stage 1
        self.beginIfBlock("count <= " + str(self.step_to_begin_adding_boost_potential))
        
        # Conventional / aMD Run
        self.addComputePerDof("sigma", "sqrt(kT/m)")
        self.addComputePerDof("fg", "f"+str(self.group))
        
        # V step
        self.addComputePerDof("v", "v + 0.5*dt*fprime/m;fprime=(f+fg*((1.0-modify) + modify*(alpha/(alpha+E-groupEnergy))^2)-fg); modify=step(E-groupEnergy)")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()

        # O step
        self.addComputePerDof("v", "a*v + b*sigma*gaussian")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()
        
        self.addComputeGlobal("groupEnergy", "energy"+str(self.group))
        self.addComputeGlobal("currentEnergy", "groupEnergy")
        
        # V step
        self.addComputePerDof("v", "v + 0.5*dt*fprime/m;fprime=(f+fg*((1.0-modify) + modify*(alpha/(alpha+E-groupEnergy))^2)-fg); modify=step(E-groupEnergy)")
        self.addConstrainVelocities()
        
        self.endBlock()
        
    def add_stage_2_and_3_instructions(self):
        super(GamdGroupBoostPotentialLangevinIntegrator, self).add_stage_2_and_3_instructions()
        #
        # Do GaMD Steps to calculate the boostPotential and the scaling Force
        # and then update the position and velocity
        #
        self.addUpdateContextState()
        
        self.addComputeGlobal("boostPotential", "0.5 * k0 * (E-groupEnergy)^2/(Vmax-Vmin)")
        self.addComputeGlobal("groupForceScalingFactor", "1.0-((k0 * (E - groupEnergy))/(Vmax - Vmin)) ")
        self.addComputePerDof("fg", "f"+str(self.group))
        
        self.addComputePerDof("sigma", "sqrt(kT/m)")
        
        # V step
        self.addComputePerDof("v", "v + 0.5*(dt*(f+(fg*groupForceScalingFactor-fg))/m)")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()

        # O step
        self.addComputePerDof("v", "a*v + b*sigma*gaussian")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()
        
        self.addComputeGlobal("groupEnergy", "energy"+str(self.group))
        self.addComputeGlobal("currentEnergy", "groupEnergy")
        self.addComputeGlobal("Vmax", "max(Vmax,groupEnergy)")
        self.addComputeGlobal("Vmin", "min(Vmin,groupEnergy)")
        self.addComputeGlobal("groupForceScalingFactor", "1.0-((k0 * (E - groupEnergy))/(Vmax - Vmin)) ")
        
        # V step
        self.addComputePerDof("v", "v + 0.5*(dt*(f+(fg*groupForceScalingFactor-fg))/m)")
        self.addConstrainVelocities()
        
        
        self.endBlock()  # count >= step_to_begin_adding_boost_potential
        
    def get_group_force_scaling_factor(self):
        return self.getGlobalVariableByName("groupForceScalingFactor")
        
    def get_group_energy(self):
        return self.getGlobalVariableByName("groupEnergy")

class GamdDualBoostPotentialLangevinIntegrator(GamdIntegratorBase):

    def __init__(self, dt, ntcmd, nteb, ntave, sigma0, sigma0Group, temperature, gamma, group):
        self.group = group
        self.sigma0Group = sigma0Group
        super(GamdDualBoostPotentialLangevinIntegrator, self).__init__(dt, ntcmd, nteb, ntave, sigma0)
        kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA
        a = numpy.exp(-gamma*dt) # current velocity component
        b = numpy.sqrt(1 - numpy.exp(- 2 * gamma * dt)) # random velocity component
        #b = numpy.sqrt(kt*(1-gamma*gamma))
        self.kT = kB * temperature
        self.a = a
        self.b = b
        self.gamma = gamma
        self.addGlobalVariable("kT", self.kT) # thermal energy
        self.addGlobalVariable("a", self.a)
        self.addGlobalVariable("b", self.b)
        self.addGlobalVariable("gamma", self.gamma)
        self.addPerDofVariable("sigma", 0) 
        
        
        # currentEnergy is the total energy
        self.addComputeGlobal("totalEnergy", "energy")
        self.addComputeGlobal("groupEnergy", "energy"+str(group))
        self.addComputeGlobal("currentEnergy", "energy")
        self.addComputeGlobal("Vmax", "max(Vmax,totalEnergy)")
        self.addComputeGlobal("Vmin", "min(Vmin,totalEnergy)")
        self.addComputeGlobal("VmaxGroup", "max(VmaxGroup,groupEnergy)")
        self.addComputeGlobal("VminGroup", "min(VminGroup,groupEnergy)")
        
        self.addUpdateContextState()
        self.add_stage_1_conventional_md_instructions()
        self.add_common_stage_1_and_2_instructions()
        self.add_stage_2_and_3_instructions()
        
        self.addComputePerDof("coordinates","x")
        
    def add_common_variables(self):
        super(GamdDualBoostPotentialLangevinIntegrator, self).add_common_variables()
        self.addGlobalVariable("sigma0Group", self.sigma0Group)
        self.addGlobalVariable("EGroup", -1E99)
        self.addGlobalVariable("k0Group", 0)

        self.addGlobalVariable("VmaxGroup", -1E99)
        self.addGlobalVariable("VminGroup", 1E99)
        self.addGlobalVariable("VavgGroup", 0)
        self.addGlobalVariable("oldVavgGroup", 0)
        self.addGlobalVariable("sigmaVGroup", 0)

        self.addGlobalVariable("M2Group", 0)
        self.addGlobalVariable("wVmaxGroup", -1E99)
        self.addGlobalVariable("wVminGroup", 1E99)
        self.addGlobalVariable("wVavgGroup", 0)
        self.addGlobalVariable("wVarianceGroup", 0)
        
        self.addGlobalVariable("groupForceScalingFactor", 0)
        self.addGlobalVariable("totalForceScalingFactor", 0)
        self.addGlobalVariable("groupEnergy", 0)
        #self.addGlobalVariable("totalEnergy", 0)
        self.addGlobalVariable("factorGroup", 0)
        self.addGlobalVariable("factorTotal", 0)
        self.addPerDofVariable("fg", 0)
        
    def add_stage_1_conventional_md_instructions(self):
        # Stage 1
        self.beginIfBlock("count <= " + str(self.step_to_begin_adding_boost_potential))
        
        # Conventional / aMD Run
        self.addComputePerDof("sigma", "sqrt(kT/m)")
        self.addComputePerDof("fg", "f"+str(self.group))
        self.addComputeGlobal("factorGroup", "((1.0-modify) + modify*(alpha/(alpha+EGroup-groupEnergy))^2); modify=step(EGroup-groupEnergy)")
        self.addComputeGlobal("factorTotal", "((1.0-modify) + modify*(alpha/(alpha+E-totalEnergy))^2); modify=step(E-totalEnergy)")
        
        # V step
        self.addComputePerDof("v", "v + 0.5*dt*fprime/m; fprime=(f+fg*factorGroup-fg)*factorTotal")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()

        # O step
        self.addComputePerDof("v", "a*v + b*sigma*gaussian")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()
        
        self.addComputeGlobal("groupEnergy", "energy"+str(self.group))
        self.addComputeGlobal("totalEnergy", "energy")
        
        # V step
        #self.addComputeGlobal("factor", "((1.0-modify) + modify*(alpha/(alpha+E-groupEnergy-totalEnergy))^2)")
        self.addComputeGlobal("factorGroup", "((1.0-modify) + modify*(alpha/(alpha+EGroup-groupEnergy))^2); modify=step(EGroup-groupEnergy)")
        self.addComputeGlobal("factorTotal", "((1.0-modify) + modify*(alpha/(alpha+E-totalEnergy))^2); modify=step(E-totalEnergy)")
        
        
        self.addComputePerDof("v", "v + 0.5*dt*fprime/m; fprime=(f+fg*factorGroup-fg)*factorTotal")
        self.addConstrainVelocities()
        
        self.endBlock()
        
    def add_stage_2_and_3_instructions(self):
        super(GamdDualBoostPotentialLangevinIntegrator, self).add_stage_2_and_3_instructions()
        #
        # Do GaMD Steps to calculate the boostPotential and the scaling Force
        # and then update the position and velocity
        #
        self.addUpdateContextState()
        
        self.addComputePerDof("fg", "f"+str(self.group))
        self.addComputeGlobal("groupForceScalingFactor", "1.0-((k0Group * (EGroup - groupEnergy))/(VmaxGroup - VminGroup)) ")
        self.addComputeGlobal("totalForceScalingFactor", "1.0-((k0 * (E - totalEnergy))/(Vmax - Vmin)) ")
        
        self.addComputePerDof("sigma", "sqrt(kT/m)")
        
        # V step
        self.addComputePerDof("v", "v + 0.5*(dt*(f+(fg*groupForceScalingFactor-fg))*totalForceScalingFactor/m)")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()

        # O step
        self.addComputePerDof("v", "a*v + b*sigma*gaussian")
        self.addConstrainVelocities()

        # R step
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()
        
        self.addComputeGlobal("groupEnergy", "energy"+str(self.group))
        self.addComputeGlobal("totalEnergy", "energy")
        self.addComputeGlobal("Vmax", "max(Vmax,totalEnergy)")
        self.addComputeGlobal("Vmin", "min(Vmin,totalEnergy)")
        self.addComputeGlobal("VmaxGroup", "max(VmaxGroup,groupEnergy)")
        self.addComputeGlobal("VminGroup", "min(VminGroup,groupEnergy)")
        self.addComputeGlobal("boostPotential", "0.5 * k0 * (E-totalEnergy)^2/(Vmax-Vmin) + 0.5 * k0Group * (EGroup-groupEnergy)^2/(VmaxGroup-VminGroup)")
        self.addComputeGlobal("groupForceScalingFactor", "1.0-((k0Group * (EGroup - groupEnergy))/(VmaxGroup - VminGroup)) ")
        self.addComputeGlobal("totalForceScalingFactor", "1.0-((k0 * (E - totalEnergy))/(Vmax - Vmin)) ")
        
        # V step
        self.addComputePerDof("v", "v + 0.5*(dt*(f+(fg*groupForceScalingFactor-fg))*totalForceScalingFactor/m)")
        self.addConstrainVelocities()
        
        
        self.endBlock()  # count >= step_to_begin_adding_boost_potential
        
    def get_group_force_scaling_factor(self):
        return self.getGlobalVariableByName("groupForceScalingFactor")
        
    def get_total_force_scaling_factor(self):
        return self.getGlobalVariableByName("totalForceScalingFactor")
        
    def get_group_energy(self):
        return self.getGlobalVariableByName("groupEnergy")
        
    def get_total_energy(self):
        return self.getGlobalVariableByName("totalEnergy")
        
    def update_potential_state_values_with_window_potential_state_values(self):
        super(GamdDualBoostPotentialLangevinIntegrator, self).update_potential_state_values_with_window_potential_state_values()
        # Update window variables
        self.addComputeGlobal("VavgGroup", "wVavgGroup")
        self.addComputeGlobal("sigmaVGroup", "sqrt(wVarianceGroup)")

        # reset variables
        self.addComputeGlobal("M2Group", "0")
        self.addComputeGlobal("wVmaxGroup", "-1E99")
        self.addComputeGlobal("wVminGroup", "1E99")
        self.addComputeGlobal("wVavgGroup", "0")
        self.addComputeGlobal("oldVavgGroup", "0")
        self.addComputeGlobal("wVarianceGroup", "0")
        
    def add_energy_value_calculations(self):
        super(GamdDualBoostPotentialLangevinIntegrator, self).add_energy_value_calculations()
        self.addComputeGlobal("wVmaxGroup", "max(groupEnergy, wVmaxGroup)")
        self.addComputeGlobal("wVminGroup", "min(groupEnergy, wVminGroup)")

        self.addComputeGlobal("oldVavgGroup", " wVavgGroup")
        self.addComputeGlobal("wVavgGroup", "wVavgGroup + (groupEnergy-wVavgGroup)/wcount")
        self.addComputeGlobal("M2Group", "M2Group + (groupEnergy - oldVavgGroup)*(groupEnergy - wVavgGroup)")
        self.addComputeGlobal("wVarianceGroup", "select(wcount-1,M2Group/(wcount - 1),0)")

    def add_global_variables(self):
        super(GamdDualBoostPotentialLangevinIntegrator, self).add_global_variables()
        self.addGlobalVariable("k0primeGroup", 0)
        self.addGlobalVariable("k0doubleprimeGroup", 0)

class GamdBoostPotentialIntegratorLowerBound(GamdIntegratorBase):

    def add_calculate_E_k0(self):
        self.addComputeGlobal("E", "Vmax")
        self.addComputeGlobal("k0prime", "(sigma0/sigmaV) * (Vmax - Vmin)/(Vmax - Vavg)")
        self.addComputeGlobal("k0", "min(1.0, k0prime)")


class GamdBoostPotentialIntegratorUpperBound(GamdIntegratorBase):

    def add_calculate_E_k0(self):
        self.addComputeGlobal("k0", "1.0")  # Our Else Case for the If Block
        self.addComputeGlobal("k0doubleprime", "1-(sigma0/sigmaV) * (Vmax - Vmin)/(Vavg - Vmin)")
        self.beginIfBlock("0.0 < k0doubleprime")
        self.beginIfBlock("k0doubleprime <= 1.0")
        self.addComputeGlobal("k0", "k0doubleprime")
        self.addComputeGlobal("E", "Vmin + (Vmax - Vmin)/k0")
        self.endBlock()
        self.endBlock()
        
class GamdDualBoostPotentialIntegratorLowerBound(GamdIntegratorBase):

    def add_calculate_E_k0(self):
        self.addComputeGlobal("E", "Vmax")
        self.addComputeGlobal("k0prime", "(sigma0/sigmaV) * (Vmax - Vmin)/(Vmax - Vavg)")
        self.addComputeGlobal("k0", "min(1.0, k0prime)")
        
        self.addComputeGlobal("EGroup", "VmaxGroup")
        self.addComputeGlobal("k0primeGroup", "(sigma0Group/sigmaVGroup) * (VmaxGroup - VminGroup)/(VmaxGroup - VavgGroup)")
        self.addComputeGlobal("k0Group", "min(1.0, k0primeGroup)")


class GamdDualBoostPotentialIntegratorUpperBound(GamdIntegratorBase):

    def add_calculate_E_k0(self):
        self.addComputeGlobal("k0", "1.0")  # Our Else Case for the If Block
        self.addComputeGlobal("k0doubleprime", "1-(sigma0/sigmaV) * (Vmax - Vmin)/(Vavg - Vmin)")
        self.beginIfBlock("0.0 < k0doubleprime")
        self.beginIfBlock("k0doubleprime <= 1.0")
        self.addComputeGlobal("k0", "k0doubleprime")
        self.addComputeGlobal("E", "Vmin + (Vmax - Vmin)/k0")
        self.endBlock()
        self.endBlock()
        
        self.addComputeGlobal("k0Group", "1.0")  # Our Else Case for the If Block
        self.addComputeGlobal("k0doubleprimeGroup", "1-(sigma0Group/sigmaVGroup) * (VmaxGroup - VminGroup)/(VavgGroup - VminGroup)")
        self.beginIfBlock("0.0 < k0doubleprimeGroup")
        self.beginIfBlock("k0doubleprimeGroup <= 1.0")
        self.addComputeGlobal("k0Group", "k0doubleprimeGroup")
        self.addComputeGlobal("EGroup", "VminGroup + (VmaxGroup - VminGroup)/k0Group")
        self.endBlock()
        self.endBlock()

class GamdTotalBoostPotentialLangevinIntegratorLowerBound(GamdTotalBoostPotentialLangevinIntegrator,
                                                  GamdBoostPotentialIntegratorLowerBound):
    def __init__(self, dt, ntcmd, nteb, ntave, sigma0, temperature, gamma):
        GamdTotalBoostPotentialLangevinIntegrator.__init__(self, dt, ntcmd, nteb, ntave, sigma0, 
                                                   temperature, gamma)
                                                   
class GamdTotalBoostPotentialLangevinIntegratorUpperBound(GamdTotalBoostPotentialLangevinIntegrator,
                                                  GamdBoostPotentialIntegratorUpperBound):
    def __init__(self, dt, ntcmd, nteb, ntave, sigma0, temperature, gamma):
        GamdTotalBoostPotentialLangevinIntegrator.__init__(self, dt, ntcmd, nteb, ntave, sigma0, 
                                                   temperature, gamma)

class GamdGroupBoostPotentialLangevinIntegratorLowerBound(GamdGroupBoostPotentialLangevinIntegrator,
                                                  GamdBoostPotentialIntegratorLowerBound):
    def __init__(self, dt, ntcmd, nteb, ntave, sigma0, temperature, gamma, group):
        GamdGroupBoostPotentialLangevinIntegrator.__init__(self, dt, ntcmd, nteb, ntave, sigma0, 
                                                   temperature, gamma, group)
        
class GamdGroupBoostPotentialLangevinIntegratorUpperBound(GamdGroupBoostPotentialLangevinIntegrator,
                                                  GamdBoostPotentialIntegratorUpperBound):
    def __init__(self, dt, ntcmd, nteb, ntave, sigma0, temperature, gamma, group):
        GamdGroupBoostPotentialLangevinIntegrator.__init__(self, dt, ntcmd, nteb, ntave, sigma0, 
                                                   temperature, gamma, group)
                                                   
class GamdDualBoostPotentialLangevinIntegratorLowerBound(GamdDualBoostPotentialLangevinIntegrator,
                                                  GamdDualBoostPotentialIntegratorLowerBound):
    def __init__(self, dt, ntcmd, nteb, ntave, sigma0, sigma0Group, temperature, gamma, group):
        GamdDualBoostPotentialLangevinIntegrator.__init__(self, dt, ntcmd, nteb, ntave, sigma0, sigma0Group, 
                                                   temperature, gamma, group)
        
class GamdDualBoostPotentialLangevinIntegratorUpperBound(GamdDualBoostPotentialLangevinIntegrator,
                                                  GamdDualBoostPotentialIntegratorUpperBound):
    def __init__(self, dt, ntcmd, nteb, ntave, sigma0, sigma0Group, temperature, gamma, group):
        GamdDualBoostPotentialLangevinIntegrator.__init__(self, dt, ntcmd, nteb, ntave, sigma0, sigma0Group, 
                                                   temperature, gamma, group)
                                                   
