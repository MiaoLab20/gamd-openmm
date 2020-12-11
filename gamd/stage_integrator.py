"""
gamd.py: Implements the GaMD integration method.

Portions copyright (c) 2020 University of Kansas
Authors: Matthew Copeland, Yinglong Miao
Contributors: Lane Votapka

"""

from __future__ import absolute_import
from enum import Enum

__author__ = "Matthew Copeland"
__version__ = "1.0"

from simtk.openmm import CustomIntegrator
from simtk import unit as unit
from abc import ABCMeta, ABC
from abc import abstractmethod

# ================
# Boost Types
# ================


class BoostType(Enum):
    TOTAL = "Total"
    DIHEDRAL = "DiHedral"


# ============================================================================================
# base class
# ============================================================================================

class GamdStageIntegrator(CustomIntegrator):
    __metaclass__ = ABCMeta

    """
        GamdIntegrator implements the GaMD integration algorithm, all modes

        Based on the following papers:
        * J. Chem Theory Comput. 2017, 13, 9-19  - DOI: 10.1021/acs.jctc.6b00931
        * J. Chem Theory Comput. 2015, 11, 3584-3595 - DOI: 10.1021/acs.jctc.5b00436

    """

    # def __init__(self,dt,alpha,E):
    def __init__(self, dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000,
                 ntebprep=200000, nteb=1000000, nstlim=3000000, ntave=50000):

        super(GamdStageIntegrator, self).__init__(dt)

        """
        Parameters
        ----------
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param nstlim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        """

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

        if nteb < ntave or nteb % ntave != 0:
            raise ValueError("nteb must be greater than and a multiple of ntave.")

        #
        #   Conventional MD Stages:
        #   Stage 1 - conventional MD preparatory stage:  no statistics are collected (equilibration of the system)
        #   Stage 2 - conventional MD  stage: boost parameters are collected (Vmax, Vmin, Vavg, and sigmaV)
        #
        #   GaMD Stages:
        #   Stage 3 - GaMD pre-equilibration stage:  boost potential is applied, boot parameters ARE NOT updated (fixed)
        #               Boost parameters from Stage 2 is applied
        #   Stage 4 - GaMD equilibration stage:  boost potential is applied, new boost parameters ARE updated
        #               Boost parameters from Stage 2 is applied.
        #   Stage 5 - GaMD production stage:  boost potential is applied, boost parameters ARE NOT updated (fixed)*
        #               Boost parameters from Stage 4 is applied.

        self.stage_1_start = 0
        self.stage_1_end = ntcmdprep
        self.stage_2_start = ntcmdprep + 1
        self.stage_2_end = ntcmd
        self.stage_2_last_ntave_window_start = ntcmd - ntave
        self.stage_3_start = ntcmd + 1
        self.stage_3_end = ntcmd + ntebprep
        self.stage_4_start = ntcmd + ntebprep + 1
        self.stage_4_end = ntcmd + nteb
        self.stage_5_start = ntcmd + nteb + 1
        self.stage_5_end = nstlim

        self.dt = dt
        self.ntcmdprep = ntcmdprep
        self.ntcmd = ntcmd
        self.ntebprep = ntebprep
        self.nteb = nteb
        self.nstlim = nstlim
        self.ntave = ntave

        #
        # NOTE:  This value is utilized to keep track of what the internal debug count is.  It's meant for
        # attempting to gain a greater understanding of the calculations or performing debug information.  Normally,
        # this value and the associated methods will go unused.
        #
        self.debug_counter = 0

        self.addGlobalVariable("stepCount", 0)
        self.addGlobalVariable("windowCount", 0)
        self.addGlobalVariable("stage", -1)
        self.addComputeGlobal("stepCount", "stepCount+1")

        self.addGlobalVariable("stageOneIfValueIsZeroOrNegative", 0)
        self.addGlobalVariable("stageTwoIfValueIsZeroOrNegative", 0)
        self.addGlobalVariable("stageThreeIfValueIsZeroOrNegative", 0)
        self.addGlobalVariable("stageFourIfValueIsZeroOrNegative", 0)
        self.addGlobalVariable("stageFiveIfValueIsZeroOrNegative", 0)

        self._add_common_variables()

        self.addGlobalVariable("starting_energy", 0.0)
        # self._add_debug()
        # self._add_debug_at_step(1)
        # self._add_debug_at_step(2)
        self.addUpdateContextState()
        self.addComputeGlobal("starting_energy", "energy")
        
        self.addComputeGlobal("stageOneIfValueIsZeroOrNegative", "(%s-stepCount)*(%s-stepCount)" % (self.stage_1_start, self.stage_1_end))
        self.addComputeGlobal("stageTwoIfValueIsZeroOrNegative", "(%s-stepCount)*(%s-stepCount)" % (self.stage_2_start, self.stage_2_end))
        self.addComputeGlobal("stageThreeIfValueIsZeroOrNegative", "(%s-stepCount)*(%s-stepCount)" % (self.stage_3_start, self.stage_3_end))
        self.addComputeGlobal("stageFourIfValueIsZeroOrNegative", "(%s-stepCount)*(%s-stepCount)" % (self.stage_4_start, self.stage_4_end))
        self.addComputeGlobal("stageFiveIfValueIsZeroOrNegative", "(%s-stepCount)*(%s-stepCount)" % (self.stage_5_start, self.stage_5_end))

        # self._add_debug()
        # self._add_debug_at_step(1)
        # self._add_debug_at_step(2)

        self._add_stage_one_instructions()
        self._add_stage_two_instructions()
        self._add_stage_three_instructions()
        self._add_stage_four_instructions()
        self._add_stage_five_instructions()

        # self._add_debug()
        # self._add_debug_at_step(1)
        # self._add_debug_at_step(2)

    def _add_debug_at_step(self, the_step):
        self.beginIfBlock("stepCount = " + str(the_step))
        self._add_debug()
        self.endBlock()
    #
    # Debugging Methods
    #

    def _add_global_debug(self, name, value):
        '''
            This is just a helper method to prepend the debug_counter to the front of a variable and add that
            global variable.

        :param name: This should be the name of the variable.
        :param value: The value to set the variable to.
        :return: none
        '''
        debug_name = str(self.debug_counter) + "_" + name
        self.addGlobalVariable(debug_name, 0.0)
        self.addComputeGlobal(debug_name, value)

    def _save_global_debug(self, name):
        '''
            If we just need to save off the value of a variable at an associated debug step, this method
            will do it with the same name with the debug counter prepended, so that we don't have to specify
            the value.
        :param name: The name of the value to save off.
        :return: none
        '''
        self._add_global_debug(name, name)

    def _add_per_dof_debug(self, name, value):
        '''
            This is just a helper method to prepend the debug_counter to the front of a variable and add
            that PerDof variable.

        :param name: This should be the name of the variable.
        :param value: The value to set the variable to.
        :return: none
        '''
        debug_name = str(self.debug_counter) + "_" + name
        self.addPerDofVariable(debug_name, 0.0)
        self.addComputePerDof(debug_name, value)

    def _save_per_dof_debug(self, name):
        """
            If we just need to save off the value of a variable at an associated debug step, this method
            will do it with the same name with the debug counter prepended, so that we don't have to specify
            the value.
        :param name: The name of the value to save off.
        :return: none
        """
        self._add_per_dof_debug(name, name)

    def _add_debug(self):
        """
            This method from the base class will add the stage and count variables to the debug set.  In addition,
            it also handles incrementing the debug counter, which is why all inheritors should call it last
            when they implement a sub-class.

        :return: none
        """
        self._save_global_debug("windowCount")
        self._save_global_debug("stage")
        self._save_global_debug("stepCount")

        self.debug_counter += 1

    def _get_global_debug_value(self, counter, name):
        return self.getGlobalVariableByName(str(counter) + "_" + name)

    def _get_per_dof_debug_value(self, counter, name):
        return self.getPerDofVariableByName(str(counter) + "_" + name)

    def get_debug_step(self, counter):
        results = {str(counter) + "_" + "windowCount": self._get_global_debug_value(counter, "windowCount"),
                   str(counter) + "_" + "stage": self._get_global_debug_value(counter, "stage"),
                   str(counter) + "_" + "stepCount": self._get_global_debug_value(counter, "stepCount")}

        return results

    def _add_stage_one_instructions(self):
        self.beginIfBlock("stepCount <= " + str(self.stage_1_end))
        # -------------------------------
        self.addComputeGlobal("stage", "1")
        self._add_conventional_md_instructions()
        # -------------------------------
        self.endBlock()

    def _add_stage_two_instructions(self):
        #self.beginIfBlock("stepCount >= " + str(self.stage_2_start))
        #self.beginIfBlock("stepCount <= " + str(self.stage_2_end))
        self.beginIfBlock("stageTwoIfValueIsZeroOrNegative <= 0")

        # -------------------------------
        self.addComputeGlobal("stage", "2")

        # Be aware:
        # In case of the need to reorder the code, we must finish with our aMD calculations prior to modifying
        # threshold_energy, so this step needs to be before we modify the threshold_energy on our last step in stage 2.
        #
        self._add_conventional_md_instructions()

        #
        # This function calculates the values for Vmax and Vmin, which should be updated throughout stages 2 and 4 only.
        #
        self._add_instructions_to_calculate_primary_boost_statistics()


        
        #
        # Just to make sure that our windowCount is set to correct the value for this point in our simulation.
        #
        #self.beginIfBlock("stepCount = " + str(self.stage_2_last_ntave_window_start))
        #self.addComputeGlobal("windowCount", "0")
        #self.endBlock()

        #
        # We only need to calculate the Vavg and sigmaV for the last ntave window in stage 2, since otherwise,
        # the values would just be overwritten, since they aren't ever used in stage 2.
        #
        self.beginIfBlock("stepCount >= " + str(self.stage_2_last_ntave_window_start))

        #
        # This helps us keep track of where we are in the ntave window.  We utilize this variable to keep a running
        # count of the average and the variance for the window.
        #
        #self.addComputeGlobal("windowCount", "windowCount + 1")

        self.addComputeGlobal("windowCount", "(1-delta(%d-stepCount))*windowCount + 1" % self.stage_2_last_ntave_window_start)
        #
        # These calculations help us to keep track of the running ntave window Vavg and variance.
        #
        self._add_instructions_to_calculate_secondary_boost_statistics()

        self.endBlock()

        #
        # If we are on the last step of the stage, we are also on the last ntave window, so we need to set the Vavg,
        # sigmaV, threshold_energy, and k0 (effective harmonic constant) we are going to use in stage 3.
        self.beginIfBlock("stepCount = " + str(self.stage_2_end))
        #
        # This method sets the values
        #
        self._update_potential_state_values_with_window_potential_state_values()
        self._calculate_threshold_energy_and_effective_harmonic_constant()
        #
        # We set the value of windowCount here, rather than in the update step, so that we can lock down where
        # modifications are occurring to the windowCount and count to this base class only.
        #
        self.addComputeGlobal("windowCount", "0")
        self.endBlock()

        # -------------------------------
        self.endBlock()
        #self.endBlock()

    def _add_stage_three_instructions(self):
        #self.beginIfBlock("stepCount >= " + str(self.stage_3_start))
        #self.beginIfBlock("stepCount <= " + str(self.stage_3_end))
        self.beginIfBlock("stageThreeIfValueIsZeroOrNegative <= 0")
        
        # -------------------------------
        self.addComputeGlobal("stage", "3")
        #
        # We recalculate the threshold energy and the effective harmonic constant at each step in stage 3.
        # These values shouldn't change though, since Vmax, Vmin, Vavg, sigma0, and sigmaV aren't changing.
        #
        self._calculate_threshold_energy_and_effective_harmonic_constant()

        self._add_gamd_instructions()
        # -------------------------------
        self.endBlock()
        #self.endBlock()

    def _add_stage_five_instructions(self):
        #self.beginIfBlock("stepCount >= " + str(self.stage_5_start))
        #self.beginIfBlock("stepCount <= " + str(self.stage_5_end))
        self.beginIfBlock("stageFiveIfValueIsZeroOrNegative <= 0")
        # -------------------------------
        self.addComputeGlobal("stage", "5")

        self._add_gamd_instructions()
        # -------------------------------
        self.endBlock()
        #self.endBlock()

    def _add_stage_four_instructions(self):
        #
        # Set our window count to zero, since this is what it should be at the start of stage 4.
        #
        #self.beginIfBlock("stepCount = " + str(self.stage_4_start))
        #self.addComputeGlobal("windowCount", "0")
        #self.endBlock()
        

        #self.beginIfBlock("stepCount >= " + str(self.stage_4_start))
        #self.beginIfBlock("stepCount <= " + str(self.stage_4_end))
        self.beginIfBlock("stageFourIfValueIsZeroOrNegative <= 0")
        # -------------------------------
        self.addComputeGlobal("stage", "4")
        self.addComputeGlobal("windowCount", "windowCount + 1")

        #
        # Throughout Stage 4, update our Vmin and Vmax.  This will continue on from our values for Vmin and Vmax from
        # stage 2, rather than resetting these values first.
        #
        self._add_instructions_to_calculate_primary_boost_statistics()

        #
        # These calculations help us to keep track of the running ntave window Vavg and variance.
        #
        self._add_instructions_to_calculate_secondary_boost_statistics()

        #
        # If we are at the end of the ntave window, then we need to calculate our Vavg, sigmaV, and reset our window
        # values.
        #
        self.beginIfBlock("windowCount = " + str(self.ntave))
        self._update_potential_state_values_with_window_potential_state_values()
        self.addComputeGlobal("windowCount", "0")
        self.endBlock()

        #
        # We recalculate the threshold energy and the effective harmonic constant based on the possible new Vmax and
        # Vmin at each step in stage 4.
        #
        self._calculate_threshold_energy_and_effective_harmonic_constant()

        #
        # WARNING:  We may have to move this step to be done first in the stage.  Recalculations of E can sometimes
        # cause issues depending on where it occurs within the step.
        #
        self._add_gamd_instructions()

        # -------------------------------
        self.endBlock()
        #self.endBlock()

    @abstractmethod
    def _add_common_variables(self):
        raise NotImplementedError("must implement _add_common_variables")

    @abstractmethod
    def _add_conventional_md_instructions(self):
        raise NotImplementedError("must implement _add_conventional_md_instructions")

    @abstractmethod
    def _add_gamd_instructions(self):
        raise NotImplementedError("must implement _add_gamd_instructions")

    @abstractmethod
    def _add_instructions_to_calculate_primary_boost_statistics(self):
        raise NotImplementedError("must implement _add_instructions_to_calculate_primary_boost_statistics")

    @abstractmethod
    def _add_instructions_to_calculate_secondary_boost_statistics(self):
        raise NotImplementedError("must implement _add_instructions_to_calculate_secondary_boost_statistics")

    @abstractmethod
    def _update_potential_state_values_with_window_potential_state_values(self):
        raise NotImplementedError("must implement _update_potential_state_values_with_window_potential_state_values")

    @abstractmethod
    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        raise NotImplementedError("must implement _calculate_threshold_energy_and_effective_harmonic_constant")

    @abstractmethod
    def get_force_scaling_factors(self):
        raise NotImplementedError("must implement get_force_scaling_factor")

    @abstractmethod
    def get_boost_potentials(self):
        raise NotImplementedError("must implement get_boost_potential")

    def get_stage(self):
        return self.getGlobalVariableByName("stage")

    def get_step_count(self):
        return self.getGlobalVariableByName("stepCount")

    def get_window_count(self):
        return self.getGlobalVariableByName("windowCount")

    def get_total_simulation_steps(self):
        return self.nstlim

    def get_coordinates(self):
        return self.getPerDofVariableByName("coordinates")

    def create_positions_file(self, filename):
        positions = self.get_coordinates()
        with open(filename, 'w') as file:
            file.write("particle, x, y, z\n")
            for i in range(len(positions)):
                line = str(i) + ", " + str(positions[i][0]) + ", " + str(positions[i][1]) + ", " + str(positions[i][2])
                file.write(line + "\n")

