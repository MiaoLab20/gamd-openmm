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

from openmm import CustomIntegrator
import openmm.unit as unit
from abc import ABC
from abc import abstractmethod

# ================
# Boost Types
# ================


class ComputeType(Enum):
    TOTAL = "Total"
    GROUP = "Group"


class BoostType(Enum):
    TOTAL = "Total"
    DIHEDRAL = "Dihedral"
    DUAL_TOTAL_DIHEDRAL = "DualTotalDihedral"
    NON_BONDED = "NonBonded"
    DUAL_NON_BONDED_DIHEDRAL = "DualNonBondedDihedral"


class BoostMethod(Enum):
    # A single Total Boost
    TOTAL = "Total"
    # A boost on one or more groups, where the order doesn't matter.
    GROUPS = "Groups"
    #  A Boost on a group and then the Total afterwards
    DUAL_DEPENDENT_GROUP_TOTAL = "DualDependentGroupTotal"


# ============================================================================================
# base class
# ============================================================================================


class GamdStageIntegrator(CustomIntegrator, ABC):

    """
        GamdIntegrator implements the GaMD integration algorithm, all modes

        Based on the following papers:
        * J. Chem Theory Comput. 2017, 13, 9-19  - 
            DOI: 10.1021/acs.jctc.6b00931
        * J. Chem Theory Comput. 2015, 11, 3584-3595 - 
            DOI: 10.1021/acs.jctc.5b00436

    """

    # def __init__(self,dt,alpha,E):
    def __init__(self, group_dict, boost_type, boost_method,
                 dt=2.0 * unit.femtoseconds,
                 ntcmdprep=200000, ntcmd=1000000,
                 ntebprep=200000, nteb=1000000, nstlim=3000000, ntave=50000):

        super(GamdStageIntegrator, self).__init__(dt)

        self.__group_dict = group_dict
        self.__boost_type = boost_type
        self._boost_method = boost_method

        """
        Parameters
        ----------
        :param group_dict: A dictionary whose keys are group indices,
            but whose values are the names of the group. 
            Example: {1:"dihedral"}
        :param total_boost: Whether to perform a total boost on this
            system.
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for 
            system equilibration.
        :param ntcmd:     The total number of conventional MD steps
            (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps 
            (including ntebprep). (must be a multiple of ntave)
        :param nstlim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the 
            average and sigma of potential energy (corresponds to 
            a running average window size).
        """

        """
        Params:
        ntcmd: Number of initial conventional MD steps used to 
            calculate the maximum, minimum, average, and standard 
            deviation of the system potential energies. Default is 
            one million timesteps for a 2fs timestep.
        nteb: Number of steps to equilibrate the system after adding 
            the boost potential. Default is one million timesteps for 
            a 2fs timestep.
        ntave: Stride to use for averaging for the calculation of Vmax,
            Vmin, etc. during conventional MD
        """
        CustomIntegrator.__init__(self, dt)

        if ntcmd < ntave or ntcmd % ntave != 0:
            raise ValueError(
                "ntcmd must be greater than and a multiple of ntave.")

        if nteb < ntave or nteb % ntave != 0:
            raise ValueError(
                "nteb must be greater than and a multiple of ntave.")

        #
        #   Conventional MD Stages:
        #   Stage 1 - conventional MD preparatory stage:  no statistics are 
        #        collected (equilibration of the system)
        #   Stage 2 - conventional MD  stage: boost parameters are collected 
        #        (Vmax, Vmin, Vavg, and sigmaV)
        #
        #   GaMD Stages:
        #   Stage 3 - GaMD pre-equilibration stage:  boost potential is 
        #        applied, boot parameters ARE NOT updated (fixed). Boost 
        #        parameters from Stage 2 is applied
        #   Stage 4 - GaMD equilibration stage:  boost potential is applied, 
        #        new boost parameters ARE updated. Boost parameters from 
        #        Stage 2 is applied.
        #   Stage 5 - GaMD production stage:  boost potential is applied, 
        #        boost parameters ARE NOT updated (fixed)* Boost parameters 
        #        from Stage 4 is applied.

        self.stage_1_start = 0
        self.stage_1_end = ntcmdprep
        self.stage_2_start = ntcmdprep + 1
        self.stage_2_end = ntcmd
        self.stage_2_last_ntave_window_start = (ntcmd - ntave) + 1
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
        # NOTE:  This value is utilized to keep track of what the internal 
        # debug count is.  It's meant for attempting to gain a greater 
        # understanding of the calculations or performing debug information. 
        # Normally, this value and the associated methods will go unused.
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
        self._setup_energy_values()

        # self._add_debug()
        # self._add_debug_at_step(1)
        # self._add_debug_at_step(2)
        self.addUpdateContextState()

        self.addComputeGlobal(
            "stageOneIfValueIsZeroOrNegative", 
            "(%s-stepCount)*(%s-stepCount)" % (self.stage_1_start, 
                                               self.stage_1_end))
        self.addComputeGlobal(
            "stageTwoIfValueIsZeroOrNegative", 
            "(%s-stepCount)*(%s-stepCount)" % (self.stage_2_start, 
                                               self.stage_2_end))
        self.addComputeGlobal(
            "stageThreeIfValueIsZeroOrNegative", 
            "(%s-stepCount)*(%s-stepCount)" % (self.stage_3_start, 
                                               self.stage_3_end))
        self.addComputeGlobal(
            "stageFourIfValueIsZeroOrNegative", 
            "(%s-stepCount)*(%s-stepCount)" % (self.stage_4_start, 
                                               self.stage_4_end))
        self.addComputeGlobal(
            "stageFiveIfValueIsZeroOrNegative", 
            "(%s-stepCount)*(%s-stepCount)" % (self.stage_5_start, 
                                               self.stage_5_end))

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

    def _setup_energy_values(self):
        self.add_global_variables_by_name("StartingPotentialEnergy", 0.0)
        self.__add_compute_globals_by_name("StartingPotentialEnergy", "{0}", ["energy"])

    def _add_debug_at_step(self, the_step):
        self.beginIfBlock("stepCount = " + str(the_step))
        self._add_debug()
        self.endBlock()
    #
    # Debugging Methods
    #

    def _add_global_debug(self, name, value):
        """
            This is just a helper method to prepend the debug_counter
            to the front of a variable and add that global variable.

        :param name: This should be the name of the variable.
        :param value: The value to set the variable to.
        :return: none
        """
        debug_name = str(self.debug_counter) + "_" + name
        self.addGlobalVariable(debug_name, 0.0)
        self.addComputeGlobal(debug_name, value)

    def _save_global_debug(self, name):
        """
            If we just need to save off the value of a variable at an 
            associated debug step, this method will do it with the same
            name with the debug counter prepended, so that we don't 
            have to specify the value.
        :param name: The name of the value to save off.
        :return: none
        """
        self._add_global_debug(name, name)

    def _add_per_dof_debug(self, name, value):
        """
            This is just a helper method to prepend the debug_counter
            to the front of a variable and add that PerDof variable.

        :param name: This should be the name of the variable.
        :param value: The value to set the variable to.
        :return: none
        """
        debug_name = str(self.debug_counter) + "_" + name
        self.addPerDofVariable(debug_name, 0.0)
        self.addComputePerDof(debug_name, value)

    def _save_per_dof_debug(self, name):
        """
            If we just need to save off the value of a variable at an 
            associated debug step, this method will do it with the 
            same name with the debug counter prepended, so that we 
            don't have to specify the value.
        :param name: The name of the value to save off.
        :return: none
        """
        self._add_per_dof_debug(name, name)

    def _add_debug(self):
        """
            This method from the base class will add the stage and 
            count variables to the debug set.  In addition, it also 
            handles incrementing the debug counter, which is why all 
            inheritors should call it last when they implement a 
            sub-class.

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
        results = {str(counter) + "_"
                   + "windowCount": self._get_global_debug_value(
                        counter, "windowCount"),
                   str(counter)
                   + "_" + "stage": self._get_global_debug_value(
                        counter, "stage"),
                   str(counter) + "_"
                   + "stepCount": self._get_global_debug_value(
                        counter, "stepCount")}

        return results

    def _add_stage_one_instructions(self):
        self.beginIfBlock("stepCount <= " + str(self.stage_1_end))
        # -------------------------------
        self.addComputeGlobal("stage", "1")
        self._add_conventional_md_instructions()
        # -------------------------------
        self.endBlock()

    def _add_stage_two_instructions(self):
        self.beginIfBlock("stageTwoIfValueIsZeroOrNegative <= 0")

        # -------------------------------
        self.addComputeGlobal("stage", "2")

        # Be aware:
        # In case of the need to reorder the code, we must finish with our 
        # aMD calculations prior to modifying threshold_energy, so this step 
        # needs to be before we modify the threshold_energy on our last step 
        # in stage 2.
        #
        self._add_conventional_md_instructions()

        #
        # This function calculates the values for Vmax and Vmin, which should 
        # be updated throughout stages 2 and 4 only.
        #
        if self._boost_method == BoostMethod.GROUPS \
                or self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:
            self._calculate_primary_boost_statistics(ComputeType.GROUP)

        if self._boost_method == BoostMethod.TOTAL \
                or self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:
            self._calculate_primary_boost_statistics(ComputeType.TOTAL)

        #
        # We only need to calculate the Vavg and sigmaV for the last ntave 
        # window in stage 2, since otherwise,
        # the values would just be overwritten, since they aren't ever used in 
        # stage 2.
        #
        self.beginIfBlock(
            "stepCount >= " + str(self.stage_2_last_ntave_window_start))

        #
        # This helps us keep track of where we are in the ntave window.  We 
        # utilize this variable to keep a running count of the average and 
        # the variance for the window.
        #

        self.addComputeGlobal(
            "windowCount", 
            "(1-delta(%d-stepCount))*windowCount + 1"
            % self.stage_2_last_ntave_window_start)
        #
        # These calculations help us to keep track of the running ntave window 
        # Vavg and variance.
        #
        if self._boost_method == BoostMethod.GROUPS \
                or self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:
            self._calculate_secondary_boost_statistics(ComputeType.GROUP)

        if self._boost_method == BoostMethod.TOTAL \
                or self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:
            self._calculate_secondary_boost_statistics(ComputeType.TOTAL)

        self.endBlock()

        #
        # If we are on the last step of the stage, we are also on the last 
        # ntave window, so we need to set the Vavg, sigmaV, threshold_energy, 
        # and k0 (effective harmonic constant) we are going to use in stage 3.
        self.beginIfBlock("stepCount = " + str(self.stage_2_end))
        #
        # This method sets the values
        #
        if self._boost_method == BoostMethod.GROUPS \
                or self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:

            self._update_potential_state_values_with_window_potential_state_values(ComputeType.GROUP)
            self._calculate_threshold_energy_and_effective_harmonic_constant(ComputeType.GROUP)

        if self._boost_method == BoostMethod.TOTAL \
                or self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:

            self._update_potential_state_values_with_window_potential_state_values(ComputeType.TOTAL)
            self._calculate_threshold_energy_and_effective_harmonic_constant(ComputeType.TOTAL)
        #
        # We set the value of windowCount here, rather than in the update step,
        # so that we can lock down where modifications are occurring to the 
        # windowCount and count to this base class only.
        #
        self.addComputeGlobal("windowCount", "0")
        self.endBlock()

        # -------------------------------
        self.endBlock()

    def _add_stage_three_instructions(self):
        self.beginIfBlock("stageThreeIfValueIsZeroOrNegative <= 0")
        
        # -------------------------------
        self.addComputeGlobal("stage", "3")
        self._do_boost_updates()
        # -------------------------------
        self.endBlock()

    def _add_stage_five_instructions(self):
        # self.beginIfBlock("stepCount >= " + str(self.stage_5_start))
        # self.beginIfBlock("stepCount <= " + str(self.stage_5_end))
        self.beginIfBlock("stageFiveIfValueIsZeroOrNegative <= 0")
        # -------------------------------
        self.addComputeGlobal("stage", "5")
        self._do_boost_updates()
        # -------------------------------
        self.endBlock()
        # self.endBlock()

    def _do_boost_updates(self):
        if self._boost_method == BoostMethod.GROUPS or \
                self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:

            self._add_gamd_pre_calc_step(ComputeType.GROUP)
            self._add_gamd_boost_calculations_step(ComputeType.GROUP)

        if self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:
            self._add_dihedral_boost_to_total_energy()

        if self._boost_method == BoostMethod.TOTAL or \
                self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:

            self._add_gamd_pre_calc_step(ComputeType.TOTAL)
            self._add_gamd_boost_calculations_step(ComputeType.TOTAL)

        self._add_gamd_update_step()
        self._add_boosts_to_starting_energies()

    def _add_stage_four_instructions(self):

        self.beginIfBlock("stageFourIfValueIsZeroOrNegative <= 0")
        # -------------------------------
        self.addComputeGlobal("stage", "4")
        self.addComputeGlobal("windowCount", "windowCount + 1")

        #
        # NOTE:  The order of these instructions is important.
        #
        self._do_boost_updates()
        self._stage_4_boost_parameters_updates()
        # -------------------------------

        self.beginIfBlock("windowCount = " + str(self.ntave))
        self.addComputeGlobal("windowCount", "0")
        self.endBlock()

        self.endBlock()

    def _stage_4_boost_parameters_updates(self):

        if self._boost_method == BoostMethod.GROUPS or \
                self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:

            #
            # Throughout Stage 4, update our Vmin and Vmax.  This will continue
            # on from our values for Vmin and Vmax from stage 2, rather than
            # resetting these values first.
            #
            self._calculate_primary_boost_statistics(ComputeType.GROUP)

            #
            # These calculations help us to keep track of the running ntave window
            # Vavg and variance.
            #
            self._calculate_secondary_boost_statistics(ComputeType.GROUP)

            #
            # If we are at the end of the ntave window, then we need to calculate
            # our Vavg, sigmaV, and reset our window
            # values.
            #
            self.beginIfBlock("windowCount = " + str(self.ntave))
            self._update_potential_state_values_with_window_potential_state_values(
                ComputeType.GROUP)
            self.endBlock()

            #
            # We recalculate the threshold energy and the effective harmonic
            # constant based on the possible new Vmax and
            # Vmin at each step in stage 4.
            #
            self._calculate_threshold_energy_and_effective_harmonic_constant(ComputeType.GROUP)

            #
            # WARNING:  We may have to move this step to be done first in the
            # stage.  Recalculations of E can sometimes cause issues depending on
            # where it occurs within the step.
            #

        if self._boost_method == BoostMethod.TOTAL or \
                self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:

            self._calculate_primary_boost_statistics(ComputeType.TOTAL)
            self._calculate_secondary_boost_statistics(ComputeType.TOTAL)
            self.beginIfBlock("windowCount = " + str(self.ntave))
            self._update_potential_state_values_with_window_potential_state_values(
                ComputeType.TOTAL)
            self.endBlock()
            self._calculate_threshold_energy_and_effective_harmonic_constant(ComputeType.TOTAL)

    @abstractmethod
    def get_effective_harmonic_constants(self):
        raise NotImplementedError("must implement get_effective_harmonic_constants")

    @abstractmethod
    def _add_common_variables(self):
        raise NotImplementedError("must implement _add_common_variables")

    @abstractmethod
    def _add_conventional_md_instructions(self):
        raise NotImplementedError(
            "must implement _add_conventional_md_instructions")

    @abstractmethod
    def _calculate_primary_boost_statistics(self, compute_type):
        raise NotImplementedError(
            "must implement _calculate_primary_boost_statistics")

    @abstractmethod
    def _calculate_secondary_boost_statistics(self,
                                              compute_type):
        raise NotImplementedError(
            "must implement _calculate_secondary_boost_statistics")

    @abstractmethod
    def _update_potential_state_values_with_window_potential_state_values(self,
                                                                          compute_type):
        raise NotImplementedError(
            "must implement " +
            "_update_potential_state_values_with_window_potential_state_values")

    @abstractmethod
    def _calculate_threshold_energy_and_effective_harmonic_constant(self,
                                                                    compute_type):
        raise NotImplementedError(
            "must implement " +
            "_calculate_threshold_energy_and_effective_harmonic_constant")

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
                line = str(i) + ", " + str(positions[i][0]) + ", " \
                    + str(positions[i][1]) + ", " + str(positions[i][2])
                file.write(line + "\n")

    #
    # The following methods are our utility methods for managing which kind of 
    # boost we are trying to perform.
    #
    #
    #

    # This method will append a unique group name to the end of the variable.
    #
    @staticmethod
    def _append_group_name(name, group_name):
        return name + "_" + str(group_name)

    # This method will append a unique group name to the end of the variable 
    # based on the type specified.
    #
    def _append_group_name_by_type(self, name, boost_type):
        return str(name + "_" + self._get_group_name_by_type(boost_type))

    # This method will append the group variable to the string. It is primarily
    # used for referencing system names. We use _append_group_name for 
    # referencing values we are creating.
    #
    @staticmethod
    def _append_group(name, group_id):
        return name + str(group_id)

    @staticmethod
    def _get_group_name_by_type(boost_type):
        return str(boost_type.value)

    def get_variable_name_by_type(self, boost_type, name):
        return self._append_group_name_by_type(name, boost_type)
    
    def add_global_variables_by_name(self, name, value):
        for group_id in self.__group_dict:
            group_name = self.__group_dict[group_id]
            var_name = self._append_group_name(name, group_name)
            self.addGlobalVariable(var_name, value)
#            print("Registered Variable ", var_name, ": ", value)

        if self._boost_method == BoostMethod.TOTAL or \
                self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:
            var_name = self._append_group_name(name, BoostType.TOTAL.value)
            self.addGlobalVariable(var_name, value)
#            print("Registered Variable ", var_name, ": ", value)
            
        return

    #
    # Add Compute Methods
    #

    def __add_compute_global_total(self, name, expression, format_list,
                                   value_by_number=False):
        var_name = self._append_group_name(name, BoostType.TOTAL.value)
        if value_by_number:
            formatted_expression = expression.format(*format_list)
        else:
            new_formats = [
                self._append_group_name(var, BoostType.TOTAL.value)
                for var in format_list]
            formatted_expression = expression.format(*new_formats)
        self.addComputeGlobal(var_name, formatted_expression)

        return

    def __add_compute_global_group(self, name, expression, format_list,
                                   value_by_number=False):
        for group_id in self.__group_dict:
            group_name = self.__group_dict[group_id]
            var_name = self._append_group_name(name, group_name)
            if value_by_number:
                new_formats = [self._append_group(var, group_id)
                               for var in format_list]
                formatted_expression = expression.format(*new_formats)
            else:
                new_formats = [self._append_group_name(var, group_name)
                               for var in format_list]
                formatted_expression = expression.format(*new_formats)
            self.addComputeGlobal(var_name, formatted_expression)

        return

    def set_global_by_name_to_value(self, name, value, compute_type):
        expression = value
        format_list = []
        value_by_number = False
        """
            This method will allow you to specify the compute type for which
            you want to run the calculations, using the group list provided
            during instantiation.
        """

        if compute_type == ComputeType.GROUP:
            self.__add_compute_global_group(name, expression, format_list,
                                            value_by_number)
        if compute_type == ComputeType.TOTAL:
            self.__add_compute_global_total(name, expression, format_list,
                                            value_by_number)
        return

    def add_compute_global_by_name(self, name, expression, format_list,
                                   compute_type):
        """
            This method will allow you to specify the compute type for which
            you want to run the calculations, using the group list provided
            during instantiation.
        """
        # Since we never set this to true, except when we were setting a value,
        # I've moved out of the parameter list, since a new function took
        # it's place.
        value_by_number = False

        if compute_type == ComputeType.GROUP:
            self.__add_compute_global_group(name, expression, format_list,
                                            value_by_number)
        if compute_type == ComputeType.TOTAL:
            self.__add_compute_global_total(name, expression, format_list,
                                            value_by_number)
        return

    def __add_compute_globals_by_name(self, name, expression, format_list,
                                      value_by_number=True):
        """
            This method will compute all of the globals based on the
            boost method assigned during instantiation.

            As this tends to be used to set default variables to be equal
            to OpenMM values, the value_by_number has been defaulted to True
            here.

        """
        if self._boost_method == BoostMethod.TOTAL or \
                self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:

            self.__add_compute_global_total(name, expression, format_list,
                                            value_by_number)

        if self._boost_method == BoostMethod.GROUPS or \
                self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:

            self.__add_compute_global_group(name, expression, format_list,
                                            value_by_number)

    def get_group_dict(self):
        return self.__group_dict

    def _add_dihedral_boost_to_total_energy(self):
        total_energy_name = self._append_group_name("StartingPotentialEnergy",
                                                    BoostType.TOTAL.value)

        for group_id in self.__group_dict:
            group_name = self.__group_dict[group_id]
            group_boost_name = self._append_group_name("BoostPotential",
                                                       group_name)
            expression = "{0} + {1}".format(total_energy_name, group_boost_name)

            self.addComputeGlobal(total_energy_name, expression)

        return

    def _add_boosts_to_starting_energies(self):
        if (self._boost_method == BoostMethod.TOTAL or
                self._boost_method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL):

            total_energy_name = self._append_group_name("StartingPotentialEnergy",
                                                    BoostType.TOTAL.value)
            total_boost_name = self._append_group_name("BoostPotential",
                                                    BoostType.TOTAL.value)
            expression = "{0} + {1}".format(total_energy_name, total_boost_name)
            self.addComputeGlobal(total_energy_name, expression)

        return
