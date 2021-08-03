"""
gamd.py: Implements the GaMD integration method.

Portions copyright (c) 2020 University of Kansas
Authors: Matthew Copeland, Yinglong Miao
Contributors: Lane Votapka

"""

from __future__ import absolute_import

__author__ = "Matthew Copeland"
__version__ = "1.0"

from simtk import unit as unit
from abc import ABC
from abc import abstractmethod
from ..stage_integrator import GamdStageIntegrator
from ..stage_integrator import BoostType


class GamdLangevinIntegrator(GamdStageIntegrator, ABC):

    def __init__(self, group_dict, total_boost,
                 dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000,
                 ntebprep=200000, nteb=1000000, nstlim=3000000, ntave=50000,
                 collision_rate=1.0 / unit.picoseconds,
                 temperature=298.15 * unit.kelvin,
                 restart_filename=None):
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
             (including ntcmdprep). (must be multiple of ntave)
         :param ntebprep:  The number of GaMD pre-equilibration steps.
         :param nteb:      The number of GaMD equilibration steps 
             (including ntebprep). (must be a multiple of ntave)
         :param nstlim:    The total number of simulation steps.
         :param ntave:     The number of steps used to smooth the
             average and sigma of potential energy (corresponds to a 
             running average window size).
         :param collision_rate:      Collision rate (gamma) compatible
             with 1/picoseconds, default: 1.0/unit.picoseconds
         :param temperature:         "Bath" temperature value 
             compatible with units.kelvin, default: 298.15*unit.kelvin
         :param restart_filename:    The file name of the restart file.
             (default=None indicates new simulation.)
         """

        self.collision_rate = collision_rate  # gamma
        self.temperature = temperature
        self.restart_filename = restart_filename
        self.kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
        self.thermal_energy = self.kB * self.temperature  # kT

        #
        # Generally, I'm trying to put variables here that I know will be used
        # across all implementations WITHOUT the name being overloaded to have
        # another meaning for an object that inherits from this base class.  
        # No guarantee I got it perfectly correct, but that is the idea.
        #
        self.global_variables = {
            "thermal_energy": self.thermal_energy,
            "collision_rate": self.collision_rate,
            "vscale": 0.0, "fscale": 0.0,
            "noisescale": 0.0
            }

        self.per_dof_variables = {"sigma": 0}

        #
        # We need to run our super classes constructor last, since it's going 
        # to execute our other methods, which have dependencies on our 
        # variables above being setup.
        #
        super(GamdLangevinIntegrator, self).__init__(
            group_dict, total_boost, dt, ntcmdprep, ntcmd, ntebprep, nteb,
            nstlim, ntave)

    def _add_common_variables(self):
        garbage = {self.addGlobalVariable(key, value) \
                   for key, value in self.global_variables.items()}
        garbage = {self.addPerDofVariable(key, value) \
                   for key, value in self.per_dof_variables.items()}
    
    @abstractmethod
    def _add_conventional_md_pre_calc_step(self):  # O Step
        raise NotImplementedError(
            "must implement _add_conventional_md_pre_calc_step")

    @abstractmethod
    def _add_conventional_md_update_step(self):
        raise NotImplementedError(
            "must implement _add_conventional_md_update_step")

    @abstractmethod
    def _add_gamd_update_step(self):
        raise NotImplementedError("must implement _add_gamd_update_step")
    
    @abstractmethod
    def _add_gamd_pre_calc_step(self):
        raise NotImplementedError("must implement _add_gamd_pre_calc_step")

    @abstractmethod
    def _add_gamd_boost_calculations_step(self):
        raise NotImplementedError(
            "must implement _add_gamd_boost_calculations_step")

    @abstractmethod
    def _add_instructions_to_calculate_primary_boost_statistics(self):
        raise NotImplementedError(
            "must implement "\
            "_add_instructions_to_calculate_primary_boost_statistics")

    @abstractmethod
    def _add_instructions_to_calculate_secondary_boost_statistics(self):
        raise NotImplementedError(
            "must implement "\
            "_add_instructions_to_calculate_secondary_boost_statistics")

    def _add_conventional_md_instructions(self):
        self._add_conventional_md_pre_calc_step()
        self._add_conventional_md_update_step()

    def _add_gamd_instructions(self, group_only=False, total_only=False):
        self._add_gamd_pre_calc_step(group_only=group_only, total_only=total_only)
        self._add_gamd_boost_calculations_step(group_only=group_only, total_only=total_only)
        self._add_gamd_update_step(group_only=group_only, total_only=total_only)
    #
    # Debugging Methods
    #

    @staticmethod
    def _get_debug_values_as_dictionary(
            dictionary, counter, function_to_retrieve_value):
        results = {}
        for key, value in dictionary.items():
            results[str(counter) + "_" + key] = \
                function_to_retrieve_value(counter, key)
        return results

    def _add_debug(self):
        garbage = {self._save_global_debug(key) \
                   for key, value in self.global_variables.items()}
        garbage = {self._save_per_dof_debug(key) \
                   for key, value in self.per_dof_variables.items()}

        super(GamdLangevinIntegrator, self)._add_debug()

    def get_debug_step(self, counter):
        results = super(GamdLangevinIntegrator, self).get_debug_step(counter)
        results.update(self._get_debug_values_as_dictionary(
            self.global_variables, counter, self._get_global_debug_value))
        results.update(self._get_debug_values_as_dictionary(
            self.per_dof_variables, counter, self._get_per_dof_debug_value))
        return results

#
#  This integrator is the basis for all of our single boost type integrators
#  to perform them in a generic way that will work across boost types.
#


class GroupBoostIntegrator(GamdLangevinIntegrator, ABC):
    """ 
    This class is an OpenMM Integrator for doing the dihedral boost 
    for Gaussian accelerated molecular dynamics.
    """

    def __init__(self, group_dict, total_boost, dt, ntcmdprep, ntcmd, ntebprep,
                 nteb, nstlim, ntave, sigma0p, sigma0D, collision_rate, temperature, 
                 restart_filename):
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
            average and sigma of potential energy (corresponds to a 
            running average window size).
        :param sigma0p:    The upper limit of the standard deviation of the 
            potential boost that allows for accurate reweighting. Total boost
            portion.
        :param sigma0D:    The upper limit of the standard deviation of the 
            potential boost that allows for accurate reweighting. Dihedral
            boost portion.
        :param collision_rate:      Collision rate (gamma) compatible
            with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible
            with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.
            (default=None indicates new simulation.)
        """
        #
        # These variables are generated per type of boost being performed
        #
        self.global_variables_by_boost_type = {
            "Vmax": -1E99, "Vmin": 1E99,  "Vavg": 0,
            "oldVavg": 0, "sigmaV": 0, "M2": 0, "wVavg": 0, "k0": 0,
            "k0prime": 0, "k0doubleprime": 0, "k0doubleprime_window": 0,
            "boosted_energy": 0, "check_boost": 0,
            "threshold_energy": -1E99}
        #
        # These variables are always kept for reporting, regardless of boost 
        # type
        #

        self.boost_global_variables = {}
        self.boost_per_dof_variables = {"newx": 0, "coordinates": 0}
        self.debug_per_dof_variables = []

        # self.debug_per_dof_variables = ["x", "v", "f", "m"]
        self.debug_global_variables = [
            "dt", "energy", "energy0", "energy1", "energy2", "energy3", 
            "energy4"]
        self.sigma0p = sigma0p
        self.sigma0D = sigma0D
        self.debuggingIsEnabled = True

        super(GroupBoostIntegrator, self).__init__(
            group_dict, total_boost, dt, ntcmdprep, ntcmd, ntebprep, nteb, 
            nstlim, ntave, collision_rate, temperature, restart_filename)

        #
        # We have to set this value separate from the others, so that when we 
        # do a non-total boost, we will still have a total boost to report 
        # back.  In that condition, the above ForceScalingFactor will get setup
        # for appropriate boost type.
        #
        # NOTE:  THIS VALUE WILL NEED TO BE FIXED SOMEHOW FOR DUAL BOOST.
        #

        self.addGlobalVariablesByName("ForceScalingFactor", 1.0)
        self.addGlobalVariablesByName("BoostPotential", 0.0)

        # This is to make sure that we get the value for the
        # beginningPotentialEnergy at the end of each simulation step.
        # self.addGlobalVariable("beginningPotentialEnergy", 0)

        self.addComputePerDof("coordinates", "x")
        return

    #
    #
    #
    #
    # def get_starting_energy(self):
    # return self.getGlobalVariableByName("starting_energy")

    # def get_current_state(self):
    # results = {"step": self.getGlobalVariableByName("stepCount")}

    # return results
    # pass

    def _add_common_variables(self):
        unused_return_values = {self.addGlobalVariable(key, value) \
                                for key, value in
                                self.boost_global_variables.items()}
        unused_return_values = {self.addPerDofVariable(key, value) \
                                for key, value in
                                self.boost_per_dof_variables.items()}
        unused_return_values = {self.addGlobalVariablesByName(key, value) \
                                for key, value in
                                self.global_variables_by_boost_type.items()}

         #hacky?
        self.addGlobalVariable("sigma0_Total", self.sigma0p)
        self.addGlobalVariable("sigma0_Dihedral", self.sigma0D)

        super(GroupBoostIntegrator, self)._add_common_variables()
        return

    def _update_potential_state_values_with_window_potential_state_values(
            self, group_only=False, total_only=False):
        # Update window variables
        self.addComputeGlobalByName("Vavg", "{0}", ["wVavg"],
                                    group_only=group_only, 
                                    total_only=total_only)
        self.addComputeGlobalByName(
            "sigmaV", "sqrt({0}/(windowCount-1))", ["M2"],
                                    group_only=group_only, 
                                    total_only=total_only)
        
        
        # Reset variables
        self.addComputeGlobalByName("M2", "0",
                                    group_only=group_only, 
                                    total_only=total_only)
        self.addComputeGlobalByName("wVavg", "0.0",
                                    group_only=group_only, 
                                    total_only=total_only)
        self.addComputeGlobalByName("oldVavg", "0.0",
                                    group_only=group_only, 
                                    total_only=total_only)
        return

    def _add_instructions_to_calculate_primary_boost_statistics(
            self, group_only=False, total_only=False):
        self.addComputeGlobalByName("Vmax", "max({0}, {1})", 
                                    ["StartingPotentialEnergy", "Vmax"],
                                    group_only=group_only, 
                                    total_only=total_only)
        self.addComputeGlobalByName("Vmin", "min({0}, {1})", 
                                    ["StartingPotentialEnergy", "Vmin"],
                                    group_only=group_only, 
                                    total_only=total_only)
        return
        
    def _add_instructions_to_calculate_secondary_boost_statistics(
            self, group_only=False, total_only=False):
        #
        # The following calculations are used to calculate the average and 
        # variance/standard deviation, rather than calculating the average at 
        # the ntave % 0 step
        #
        # Algorithm Description:
        #
        # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance\
        #Welford's_online_algorithm
        #
        #
        
        self.addComputeGlobalByName("oldVavg", "{0}", ["wVavg"],
                                    group_only=group_only, 
                                    total_only=total_only)
        self.addComputeGlobalByName(
            "wVavg", "{0} + ({1}-{0})/windowCount", 
            ["wVavg", "StartingPotentialEnergy"],
            group_only=group_only, total_only=total_only)
        self.addComputeGlobalByName(
            "M2", "{0} + ({1}-{2})*({1}-{3})", 
            ["M2", "StartingPotentialEnergy", "oldVavg", "wVavg"],
            group_only=group_only, total_only=total_only)
        return
    
    def _add_conventional_md_pre_calc_step(self):
        self.addComputeGlobal("vscale", "exp(-dt*collision_rate)")
        self.addComputeGlobal("fscale", "(1-vscale)/collision_rate")
        self.addComputeGlobal("noisescale", 
                              "sqrt(thermal_energy*(1-vscale*vscale))")
        return

    def _add_conventional_md_update_step(self):
        self.addComputePerDof("newx", "x")
        self.addComputePerDof("v", 
                              "vscale*v + fscale*f/m "\
                              "+ noisescale*gaussian/sqrt(m)")
        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-newx)/dt")
        return

    def _add_gamd_pre_calc_step(self, group_only=False, total_only=False, 
                                skip_scales=False):
        if not skip_scales:
            self.addComputeGlobal("vscale", "exp(-dt*collision_rate)")
            self.addComputeGlobal("fscale", "(1-vscale)/collision_rate")
            self.addComputeGlobal("noisescale", 
                                  "sqrt(thermal_energy*(1-vscale*vscale))")
        #
        # We do not apply the boost potential to the energy value since 
        # energy is read only.
        #
        
        self.addComputeGlobalByName(
            "BoostPotential", "0.5 * {0} * ({1} - {2})^2 / ({3} - {4})", 
            ["k0", "threshold_energy", "StartingPotentialEnergy", "Vmax", 
             "Vmin"], group_only=group_only, total_only=total_only)
        
        
        
        #
        # "BoostPotential*step(threshold_energy-boosted_energy)")
        self.addComputeGlobalByName(
            "BoostPotential", "{0}*step({1} - ({2} + {3}))", 
            ["BoostPotential", "threshold_energy", "BoostPotential", 
             "StartingPotentialEnergy"], group_only=group_only, 
            total_only=total_only)
        
        #
        # If the boostPotential is zero, we want to set the Force Scaling 
        # Factor to one, which is what we will use the check_boost value to 
        # do in a later portion of the code.
        #
        self.addComputeGlobalByName(
            "check_boost", "1 - delta({0})", ["BoostPotential"], 
            group_only=group_only, total_only=total_only)

        # "boosted_energy" = "energy + BoostPotential"
        self.addComputeGlobalByName(
            "boosted_energy", "{0} + {1}", 
            ["StartingPotentialEnergy", "BoostPotential"], 
            group_only=group_only, total_only=total_only)
        
        return
        
    def _add_gamd_boost_calculations_step(self, group_only=False, total_only=False):
        
        self.addComputeGlobalByName(
            "ForceScalingFactor", "1.0 - (({0} * ({1} - {2}))/({3} - {4}))", 
            ["k0", "threshold_energy", "StartingPotentialEnergy", "Vmax", 
             "Vmin"], group_only=group_only, total_only=total_only)
        
        # This is the psuedo code of what we are about to do, in case it helps 
        # you read it.
        #
        # self.beginIfBlock("boosted_energy >= threshold_energy")
        #

        #
        #  When the boosted energy is greater than or equal to the threshold 
        # energy, the value of check_boost will be 0. This will cause the 
        # following equation to change the ForceScalingFactor to 1.0.  When 
        # the boosted_energy is less than the threshold energy, we are in our 
        # normal good condition, and just want to keep the ForceScalingFactor 
        # the same.
        #
        #   NOTE:  We do these odd computational gymnastics to counteract the 
        # problem within OpenMM with if statements causing the JIT compiler to 
        # take an exponentially larger amount of time to start.
        #
        #   1.0 - 1.0 * check_boost + check_boost * ForceScalingFactor"
        self.addComputeGlobalByName(
            "ForceScalingFactor", "1.0 - {0} + {0} * {1}", 
            ["check_boost", "ForceScalingFactor"], group_only=group_only, 
            total_only=total_only)
        
        #
        #
        #
        return

    def _add_gamd_update_step(self, group_only=False, total_only=False):
        self.addComputePerDof("newx", "x")
        #
        if self.get_total_boost() and not group_only:
            if len(self.get_group_dict()) == 0: # only Total boost
                self.addComputePerDof(
                    "v", "vscale*v + fscale*f*{0}/m "\
                    "+ noisescale*gaussian/sqrt(m)".format(
                        self._append_group_name(
                            "ForceScalingFactor", BoostType.TOTAL.value)))
            else:
                self.addComputePerDof(
                    "v", "vscale*v + fscale*f0*{0}/m "\
                    "+ noisescale*gaussian/sqrt(m)".format(
                        self._append_group_name(
                            "ForceScalingFactor", BoostType.TOTAL.value)))
            
            if not total_only:
                for group_id in self.get_group_dict():
                    group_name = self.get_group_dict()[group_id]
                    self.addComputePerDof("v", "v + fscale*{0}*{1}*{2}/m"
                        .format(self._append_group("f", group_id), 
                                self._append_group_name("ForceScalingFactor", 
                                                        BoostType.TOTAL.value),
                                self._append_group_name("ForceScalingFactor", 
                                                        group_name)))
        
        elif not total_only:
            self.addComputePerDof(
                "v", "vscale*v + fscale*f0/m + noisescale*gaussian/sqrt(m)")
            for group_id in self.get_group_dict():
                group_name = self.get_group_dict()[group_id]
                self.addComputePerDof("v", "v + fscale*{0}*{1}/m"
                    .format(self._append_group("f", group_id), 
                            self._append_group_name("ForceScalingFactor", 
                                                    group_name)))
        
        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-newx)/dt")
        return

    def get_force_scaling_factors(self):
        force_scaling_factors = {
            self._append_group_name(
                "ForceScalingFactor", BoostType.TOTAL.value): 1.0,
            self._append_group_name(
                "ForceScalingFactor", "Dihedral"): 1.0
                                 }
        
        for group_id in self.get_group_dict():
            group_name = self.get_group_dict()[group_id]
            var_name = self._append_group_name("ForceScalingFactor", 
                                               group_name)
            var_value = self.getGlobalVariableByName(var_name)
            force_scaling_factors[var_name] = var_value
        
        if self.get_total_boost():
            var_name = self._append_group_name("ForceScalingFactor", 
                                               BoostType.TOTAL.value)
            var_value = self.getGlobalVariableByName(var_name)
            force_scaling_factors[var_name] = var_value
        
        return force_scaling_factors

    def get_boost_potentials(self):
        boost_potentials = {
            self._append_group_name(
                "BoostPotential", BoostType.TOTAL.value): 0.0,
            self._append_group_name(
                "BoostPotential", "Dihedral"): 0.0
                                 }
        for group_id in self.get_group_dict():
            group_name = self.get_group_dict()[group_id]
            var_name = self._append_group_name("BoostPotential", 
                                               group_name)
            var_value = self.getGlobalVariableByName(var_name)
            boost_potentials[var_name] = var_value
            
        if self.get_total_boost():
            var_name = self._append_group_name("BoostPotential", 
                                               BoostType.TOTAL.value)
            var_value = self.getGlobalVariableByName(var_name)
            boost_potentials[var_name] = var_value
            
        
                        
        return boost_potentials

    def __calculate_simple_threshold_energy_and_effective_harmonic_constant_total(
            self):
        self.addComputeGlobalByName("threshold_energy", "{0}",["Vmax"], total_only=True)
        self.addComputeGlobalByName(
            "k0prime", "({0}/{1}) * ({2} - {3}) / ({2} - {4})",
            ["sigma0", "sigmaV", "Vmax", "Vmin", "Vavg"], total_only=True)
        self.addComputeGlobalByName("k0", "min(1.0, {0})", ["k0prime"], total_only=True)
        
        return
    
    def __calculate_simple_threshold_energy_and_effective_harmonic_constant_group(
            self, group_id):
        self.addComputeGlobalByName("threshold_energy", "{0}",["Vmax"], group_only=True)
        self.addComputeGlobalByName(
            "k0prime", "({0}/{1}) * ({2} - {3}) / ({2} - {4})",
            ["sigma0", "sigmaV", "Vmax", "Vmin", "Vavg"], group_only=True)
        self.addComputeGlobalByName("k0", "min(1.0, {0})", ["k0prime"], group_only=True)
        return
        
    def _upper_bound_calculate_threshold_energy_and_effective_harmonic_constant(
            self, group_only=False, total_only=False):
        self.addComputeGlobalByName("k0", "1.0", group_only=group_only, 
            total_only=total_only)
        self.addComputeGlobalByName(
            "k0doubleprime", "(1 - {0}/{1}) * ({2} - {3})/({4} - {3})",
            ["sigma0", "sigmaV", "Vmax", "Vmin", "Vavg"], group_only=group_only, 
            total_only=total_only)
        self.addComputeGlobalByName("k0", "{0}", ["k0doubleprime"], group_only=group_only, 
            total_only=total_only)
        self.addComputeGlobalByName(
            "threshold_energy", "{0} + ({1} - {0})/{2}", ["Vmin", "Vmax", "k0"],
            group_only=group_only, total_only=total_only)
        self.addComputeGlobalByName(
            "k0doubleprime_window", "(-{0}) * (1 - {0})", ["k0doubleprime"], 
            group_only=group_only, total_only=total_only)
        
        if self.get_total_boost() and not group_only:
            self.beginIfBlock(
                self._append_group_name("k0doubleprime_window", 
                                        BoostType.TOTAL.value) + " >= 0.0")
            self.__calculate_simple_threshold_energy_and_effective_harmonic_constant_total()
            self.endBlock()
        
        if not total_only:
            for group_id in self.get_group_dict():
                group_name = self.get_group_dict()[group_id]
                self.beginIfBlock(
                    self._append_group_name("k0doubleprime_window", group_name) \
                    + " >= 0.0")
                self.__calculate_simple_threshold_energy_and_effective_harmonic_constant_group(group_id)
                self.endBlock()
        return

    def _lower_bound_calculate_threshold_energy_and_effective_harmonic_constant(
            self, group_only=False, total_only=False):
        if not total_only:
            for group_id in self.get_group_dict():
                self.__calculate_simple_threshold_energy_and_effective_harmonic_constant_group(group_id)
        
        if not group_only:    
            self.__calculate_simple_threshold_energy_and_effective_harmonic_constant_total()
        return
