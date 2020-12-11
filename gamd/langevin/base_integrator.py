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
from abc import ABCMeta, ABC
from abc import abstractmethod
import numpy
from ..stage_integrator import GamdStageIntegrator
from ..stage_integrator import BoostType


class GamdLangevinIntegrator(GamdStageIntegrator, ABC):

    def __init__(self,
                 dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000,
                 ntebprep=200000, nteb=1000000, nstlim=3000000, ntave=50000,
                 collision_rate=1.0 / unit.picoseconds,
                 temperature=298.15 * unit.kelvin,
                 restart_filename=None):
        """
         Parameters
         ----------
         :param dt:        The Amount of time between each time step.
         :param ntcmdprep: The number of conventional MD steps for system equilibration.
         :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be multiple of ntave)
         :param ntebprep:  The number of GaMD pre-equilibration steps.
         :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
         :param nstlim:    The total number of simulation steps.
         :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to
                           a running average window size).
         :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
         :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
         :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
         """

        self.collision_rate = collision_rate  # gamma
        self.temperature = temperature
        self.restart_filename = restart_filename
        self.kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
        self.thermal_energy = self.kB * self.temperature  # kT
        #self.current_velocity_component = numpy.exp(-self.collision_rate * dt)  # a
        #self.random_velocity_component = numpy.sqrt(1 - numpy.exp(- 2 * self.collision_rate * dt))  # b

        #
        # Generally, I'm trying to put variables here that I know will be used across all implementations WITHOUT the
        # name being overloaded to have another meaning for an object that inherits from this base class.  No guarantee
        # I got it perfectly correct, but that is the idea.
        #
        self.global_variables = {"thermal_energy": self.thermal_energy,
                                 #"current_velocity_component": self.current_velocity_component,
                                 #"random_velocity_component": self.random_velocity_component,
                                 "collision_rate": self.collision_rate,
                                 "vscale": 0.0, "fscale": 0.0,
                                 "noisescale": 0.0
                                 }

        self.per_dof_variables = {"sigma": 0}

        #
        # We need to run our super classes constructor last, since it's going to execute our other methods, which
        # have dependencies on our variables above being setup.
        #
        super(GamdLangevinIntegrator, self).__init__(dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave)

    def _add_common_variables(self):
        garbage = {self.addGlobalVariable(key, value) for key, value in self.global_variables.items()}
        garbage = {self.addPerDofVariable(key, value) for key, value in self.per_dof_variables.items()}
    
    @abstractmethod
    def _add_conventional_md_pre_calc_step(self):  # O Step
        raise NotImplementedError("must implement _add_conventional_md_pre_calc_step")
    '''
    @abstractmethod
    def _add_conventional_md_position_update_step(self):  # R Step
        raise NotImplementedError("must implement _add_conventional_md_position_update_step")

    @abstractmethod
    def _add_conventional_md_velocity_update_step(self):  # V Step
        raise NotImplementedError("must implement _add_conventional_md_velocity_update_step")
    
    @abstractmethod
    def _add_conventional_md_stochastic_velocity_update_step(self):  # O Step
        raise NotImplementedError("must implement _add_conventional_md_stochastic_velocity_update_step")
    '''
    @abstractmethod
    def _add_conventional_md_update_step(self):
        raise NotImplementedError("must implement _add_conventional_md_update_step")
    '''
    @abstractmethod
    def _add_gamd_position_update_step(self):  # R Step
        raise NotImplementedError("must implement _add_gamd_position_update_step")

    @abstractmethod
    def _add_gamd_velocity_update_step(self):  # V Step
        raise NotImplementedError("must implement _add_gamd_velocity_update_step")
    
    @abstractmethod
    def _add_gamd_stochastic_velocity_update_step(self):  # O Step
        raise NotImplementedError("must implement _add_gamd_stochastic_velocity_update_step")
    '''
    @abstractmethod
    def _add_gamd_update_step(self):
        raise NotImplementedError("must implement _add_gamd_update_step")
    
    @abstractmethod
    def _add_gamd_pre_calc_step(self):
        raise NotImplementedError("must implement _add_gamd_pre_calc_step")

    @abstractmethod
    def _add_gamd_boost_calculations_step(self):
        raise NotImplementedError("must implement _add_gamd_boost_calculations_step")

    @abstractmethod
    def _add_instructions_to_calculate_primary_boost_statistics(self):
        raise NotImplementedError("must implement _add_instructions_to_calculate_primary_boost_statistics")

    @abstractmethod
    def _add_instructions_to_calculate_secondary_boost_statistics(self):
        raise NotImplementedError("must implement _add_instructions_to_calculate_secondary_boost_statistics")

    def _add_conventional_md_instructions(self):
        self._add_conventional_md_pre_calc_step()
        '''
        self._add_conventional_md_velocity_update_step()
        self._add_conventional_md_position_update_step()
        self._add_conventional_md_stochastic_velocity_update_step()
        self._add_conventional_md_position_update_step()
        self._add_conventional_md_velocity_update_step()
        '''
        self._add_conventional_md_update_step()

    def _add_gamd_instructions(self):
        self._add_gamd_pre_calc_step()
        self._add_gamd_boost_calculations_step()
        
        '''
        self._add_gamd_velocity_update_step()
        self._add_gamd_position_update_step()
        self._add_gamd_stochastic_velocity_update_step()
        self._add_gamd_position_update_step()

        #
        # We should only need to calculating the scaling factor once per step, since Vmax, Vmin, the threshold energy,
        # and the effective harmonic constant don't change after being set.  It's only a question if the energy changes
        # somehow during the step.
        #
        #self._add_gamd_boost_calculations_step()

        self._add_gamd_velocity_update_step()
        '''
        self._add_gamd_update_step()
    #
    # Debugging Methods
    #

    @staticmethod
    def _get_debug_values_as_dictionary(dictionary, counter, function_to_retrieve_value):
        results = {}
        for key, value in dictionary.items():
            results[str(counter) + "_" + key] = function_to_retrieve_value(counter, key)
        return results

    def _add_debug(self):
        garbage = {self._save_global_debug(key) for key, value in self.global_variables.items()}
        garbage = {self._save_per_dof_debug(key) for key, value in self.per_dof_variables.items()}

        super(GamdLangevinIntegrator, self)._add_debug()

    def get_debug_step(self, counter):
        results = super(GamdLangevinIntegrator, self).get_debug_step(counter)
        results.update(self._get_debug_values_as_dictionary(self.global_variables, counter, self._get_global_debug_value))
        results.update(self._get_debug_values_as_dictionary(self.per_dof_variables, counter, self._get_per_dof_debug_value))
        return results

#
#  This integrator is the basis for all of our single boost type integrators
#  to perform them in a generic way that will work across boost types.
#


class GroupBoostIntegrator(GamdLangevinIntegrator, ABC):
    """ This class is an OpenMM Integrator for doing the dihedral boost for Gaussian accelerated molecular dynamics.
    """

    def __init__(self, system_group, group_name, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, sigma0,
                 collision_rate, temperature, restart_filename):
        """
        Parameters
        ----------
        :param system_group:    This value indicates what value should be appended to system names (energy, force) for accessing the correct group's variable.
        :param group_name:  This variable along with the system_group is used to create a unique name for each of our variables, so that if you are composing groups for boosts, they do not overwrite.
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param nstlim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        :param sigma0:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """
        #
        # These variables are generated per type of boost being performed
        #
        self.global_variables_by_boost_type = {"Vmax": -1E99, "Vmin": 1E99,  "Vavg": 0,
                                               "oldVavg": 0, "sigmaV": 0, "M2": 0, "wVavg": 0, "wVariance": 0, "k0": 0,
                                               "k0prime": 0, "k0doubleprime": 0, "k0doubleprime_window": 0,
                                               "boosted_energy": 0, "check_boost": 0, "sigma0": sigma0,
                                               "threshold_energy": -1E99}
        #
        # These variables are always kept for reporting, regardless of boost type
        #

        self.boost_global_variables = {"beginningPotentialEnergy": 0}
        self.boost_per_dof_variables = {"newx": 0, "coordinates": 0}
        self.debug_per_dof_variables = []

        # self.debug_per_dof_variables = ["x", "v", "f", "m"]
        self.debug_global_variables = ["dt", "energy", "energy0", "energy1", "energy2", "energy3", "energy4"]
        self.sigma0 = sigma0
        self.debuggingIsEnabled = True

        self.__system_group = system_group
        self.__group_name = group_name

        super(GroupBoostIntegrator, self).__init__(dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, collision_rate,
                                                   temperature, restart_filename)


        #
        # We have to set this value separate from the others, so that when we do a non-total boost, we will still
        # have a total boost to report back.  In that condition, the above ForceScalingFactor will get setup for
        # appropriate boost type.
        #
        # NOTE:  THIS VALUE WILL NEED TO BE FIXED SOMEHOW FOR DUAL BOOST.
        #
        self.addGlobalVariable(self._append_group_name_by_type("ForceScalingFactor", BoostType.TOTAL), 1.0)
        self.addGlobalVariable(self._append_group_name_by_type("BoostPotential", BoostType.TOTAL), 0.0)


        if self.__group_name == BoostType.TOTAL or self.__group_name == BoostType.DIHEDRAL:
            self.addGlobalVariable(self._append_group_name_by_type("ForceScalingFactor", BoostType.DIHEDRAL), 1.0)
            self.addGlobalVariable(self._append_group_name_by_type("BoostPotential", BoostType.DIHEDRAL), 0.0)
        else:
            self.addGlobalVariable(self._append_group_name("ForceScalingFactor"), 1.0)
            self.addGlobalVariable(self._append_group_name("BoostPotential"), 0.0)




        # This is to make sure that we get the value for the beginningPotentialEnergy at the end of each simulation step.
        # self.addGlobalVariable("beginningPotentialEnergy", 0)
        self.addComputeGlobal("beginningPotentialEnergy", "energy")
        self.addComputePerDof("coordinates", "x")

    #
    # The following methods are our utility methods for managing which kind of boost we are trying to perform.
    #
    #
    #
    def _get_group_energy_name(self):
        return self._append_group("energy")

    def _get_group_name(self):
        return str(self.__group_name.value)

    @staticmethod
    def _get_group_name_by_type(boost_type):
        return str(boost_type.value)

    # This method will append a unique group name to the end of the variable.
    #
    def _append_group_name(self, name):
        return str(self._get_group_name() + name)

    # This method will append a unique group name to the end of the variable based on the type specified.
    #
    def _append_group_name_by_type(self, name, boost_type):
        return str(self._get_group_name_by_type(boost_type) + name)

    # This method will append the group variable to the string. It is primarily used for referencing system names. We
    # use __apend_group_name for referencing values we are creating.
    #
    def _append_group(self, name):
        return name + str(self.__system_group)

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
        unused_return_values = {self.addGlobalVariable(key, value) for key, value in
                                self.boost_global_variables.items()}
        unused_return_values = {self.addPerDofVariable(key, value) for key, value in
                                self.boost_per_dof_variables.items()}
        unused_return_values = {self.addGlobalVariable(self._append_group_name(key), value) for key, value in
                                self.global_variables_by_boost_type.items()}

        super(GroupBoostIntegrator, self)._add_common_variables()

    def _update_potential_state_values_with_window_potential_state_values(self):
        # Update window variables
        self.addComputeGlobal(self._append_group_name("Vavg"), self._append_group_name("wVavg"))
        self.addComputeGlobal(self._append_group_name("sigmaV"), "sqrt(" + self._append_group_name("wVariance") + ")")

        # Reset variables
        self.addComputeGlobal(self._append_group_name("M2"), "0")
        self.addComputeGlobal(self._append_group_name("wVavg"), self._append_group("energy"))
        self.addComputeGlobal(self._append_group_name("oldVavg"), self._append_group("energy"))
        self.addComputeGlobal(self._append_group_name("wVariance"), "0")

    def _add_instructions_to_calculate_primary_boost_statistics(self):
        self.addComputeGlobal(self._append_group_name("Vmax"), "max({0}, {1})".format(self._append_group("energy"),
                                                                                       self._append_group_name("Vmax")))
        self.addComputeGlobal(self._append_group_name("Vmin"), "min({0}, {1})".format(self._append_group("energy"),
                                                                                       self._append_group_name("Vmin")))

    def _add_instructions_to_calculate_secondary_boost_statistics(self):
        #
        # The following calculations are used to keep a running average,
        # rather than calculating the average at the ntave % 0 step
        #
        self.addComputeGlobal(self._append_group_name("oldVavg"), self._append_group_name("wVavg"))
        self.addComputeGlobal(self._append_group_name("wVavg"), "{0} + ({1}-{0})/windowCount".format(
            self._append_group_name("wVavg"), self._append_group("energy")))

        self.addComputeGlobal(self._append_group_name("M2"), "{0} + ({1}-{2})*({1}-{3})".format(
            self._append_group_name("M2"), self._append_group("energy"), self._append_group_name("oldVavg"),
            self._append_group_name("wVavg")))

        self.addComputeGlobal(self._append_group_name("wVariance"), "select(windowCount - 1,{0}/(windowCount - 1),0)"
                              .format(self._append_group_name("M2")))

    def _add_conventional_md_pre_calc_step(self):
        self.addComputeGlobal("vscale", "exp(-dt*collision_rate)")
        self.addComputeGlobal("fscale", "(1-vscale)/collision_rate")
        self.addComputeGlobal("noisescale", "sqrt(thermal_energy*(1-vscale*vscale))")

    def _add_conventional_md_update_step(self):
        self.addComputePerDof("newx", "x")
        self.addComputePerDof("v", "vscale*v + fscale*f/m + noisescale*gaussian/sqrt(m)")
        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-newx)/dt")

    def _add_gamd_pre_calc_step(self):
        self.addComputeGlobal("vscale", "exp(-dt*collision_rate)")
        self.addComputeGlobal("fscale", "(1-vscale)/collision_rate")
        self.addComputeGlobal("noisescale", "sqrt(thermal_energy*(1-vscale*vscale))")
        #
        # We do not apply the boost potential to the energy value since energy is read only.
        #
        self.addComputeGlobal(self._append_group_name("BoostPotential"), "0.5 * {0} * ({1} - {2})^2 / ({3} - {4})".
                              format(self._append_group_name("k0"), self._append_group_name("threshold_energy"),
                                     self._append_group("energy"), self._append_group_name("Vmax"),
                                     self._append_group_name("Vmin")))


        #
        # "BoostPotential*step(threshold_energy-boosted_energy)")
        self.addComputeGlobal(self._append_group_name("BoostPotential"), "{0}*step({1} - ({2} + {3}))".format(
            self._append_group_name("BoostPotential"), self._append_group_name("threshold_energy"),
            self._append_group_name("BoostPotential"), self._append_group("energy")))

        #
        # If the boostPotential is zero, we want to set the Force Scaling Factor to one, which is what we will use
        # the check_boost value to do in a later portion of the code.
        #
        self.addComputeGlobal(self._append_group_name("check_boost"), "1 - delta({0})".format(self._append_group_name("BoostPotential")))


        # "boosted_energy" = "energy + BoostPotential"
        self.addComputeGlobal(self._append_group_name("boosted_energy"), "{0} + {1}".format(
            self._append_group("energy"),
            self._append_group_name("BoostPotential")))

    def _add_gamd_boost_calculations_step(self):

        self.addComputeGlobal(self._append_group_name("ForceScalingFactor"), "1.0 - (({0} * ({1} - {2}))/({3} - {4}))"
                              .format(self._append_group_name("k0"), self._append_group_name("threshold_energy"),
                                      self._append_group("energy"), self._append_group_name("Vmax"),
                                      self._append_group_name("Vmin")))

        # This is the psuedo code of what we are about to do, in case it helps you read it.
        #
        # self.beginIfBlock("boosted_energy >= threshold_energy")
        #

        #
        #  When the boosted energy is greater than or equal to the threshold energy, the value of check_boost will be 0.
        #  This will cause the following equation to change the ForceScalingFactor to 1.0.  When the boosted_energy
        #  is less than the threshold energy, we are in our normal good condition, and just want to keep the
        #  ForceScalingFactor the same.
        #
        #   NOTE:  We do these odd computational gymnastics to counteract the problem within OpenMM with
        #          if statements causing the JIT compiler to take an exponentially larger amount of time to start.
        #
        #   1.0 - 1.0 * check_boost + check_boost * ForceScalingFactor"
        self.addComputeGlobal(self._append_group_name("ForceScalingFactor"), "1.0 - {0} + {0} * {1}"
                              .format(self._append_group_name("check_boost"),
                                      self._append_group_name("ForceScalingFactor")))

        #
        #
        #

    def _add_gamd_update_step(self):
        self.addComputePerDof("newx", "x")
        #
        if self.__group_name == BoostType.TOTAL:
            self.addComputePerDof("v", "vscale*v + fscale*{0}*{1}/m + noisescale*gaussian/sqrt(m)"
                                  .format(self._append_group("f"), self._append_group_name("ForceScalingFactor")))
        elif self.__group_name == BoostType.DIHEDRAL:
            # We take care of all of the forces that aren't the dihedral.
            self.addComputePerDof("v", "vscale*v + fscale*f0/m + noisescale*gaussian/sqrt(m)")
            # We boost the dihedral force.

            self.addComputePerDof("v", "v + fscale*{0}*{1}/m"
                                  .format(self._append_group("f"), self._append_group_name("ForceScalingFactor")))
        else:
            print("Failure in detecting boost type to determine proper boost methodolgy.")

        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-newx)/dt")

    def get_beginning_potential_energy(self):
        return self.getGlobalVariableByName("beginningPotentialEnergy")

    def get_force_scaling_factors(self):
        force_scaling_factors = {
            self._append_group_name_by_type("ForceScalingFactor", BoostType.TOTAL): self.getGlobalVariableByName(
                self._append_group_name_by_type("ForceScalingFactor", BoostType.TOTAL))}
        if self.__group_name == BoostType.TOTAL or self.__group_name == BoostType.DIHEDRAL:
            force_scaling_factors[self._append_group_name_by_type("ForceScalingFactor", BoostType.DIHEDRAL)] = \
                self.getGlobalVariableByName(self._append_group_name_by_type("ForceScalingFactor", BoostType.DIHEDRAL))
        else:
            force_scaling_factors[self._append_group_name("ForceScalingFactor")] = self.getGlobalVariableByName(
                self._append_group_name("ForceScalingFactor"))

        return force_scaling_factors

    def get_boost_potentials(self):
        boost_potentials = {
            self._append_group_name_by_type("BoostPotential", BoostType.TOTAL): self.getGlobalVariableByName(
                self._append_group_name_by_type("BoostPotential", BoostType.TOTAL))}

        if self.__group_name == BoostType.TOTAL or self.__group_name == BoostType.DIHEDRAL:
            boost_potentials[self._append_group_name_by_type("BoostPotential", BoostType.DIHEDRAL)] = \
                self.getGlobalVariableByName(self._append_group_name_by_type("BoostPotential", BoostType.DIHEDRAL))
        else:
            boost_potentials[self._append_group_name("BoostPotential")] = self.getGlobalVariableByName(
                self._append_group_name("BoostPotential"))

        return boost_potentials

    def __calculate_simple_threshold_energy_and_effective_harmonic_constant(self):
        self.addComputeGlobal(self._append_group_name("threshold_energy"), self._append_group_name("Vmax"))
        # "(sigma0/sigmaV) * (Vmax - Vmin)/(Vmax - Vavg)"
        self.addComputeGlobal(self._append_group_name("k0prime"),
                              "({0}/{1}) * ({2} - {3}) / ({2} - {4})".format(self._append_group_name("sigma0"),
                                                                             self._append_group_name("sigmaV"),
                                                                             self._append_group_name("Vmax"),
                                                                             self._append_group_name("Vmin"),
                                                                             self._append_group_name("Vavg")))
        self.addComputeGlobal(self._append_group_name("k0"),
                              "min(1.0, {0}) ".format(self._append_group_name("k0prime")))

    def _upper_bound_calculate_threshold_energy_and_effective_harmonic_constant(self):
        self.addComputeGlobal(self._append_group_name("k0"), "1.0")
        # "1 - (sigma0/sigmaV) * (Vmax - Vmin)/(Vavg - Vmin)"
        self.addComputeGlobal(self._append_group_name("k0doubleprime"),
                              "(1 - {0}/{1}) * ({2} - {3})/({4} - {3})".format(self._append_group_name("sigma0"),
                                                                               self._append_group_name("sigmaV"),
                                                                               self._append_group_name("Vmax"),
                                                                               self._append_group_name("Vmin"),
                                                                               self._append_group_name("Vavg")))
        #
        #
        #
        #
        self.addComputeGlobal(self._append_group_name("k0"), self._append_group_name("k0doubleprime"))
        # "Vmin + (Vmax - Vmin)/k0"
        self.addComputeGlobal(self._append_group_name("threshold_energy"),
                              "{0} + ({1} - {0})/{2}".format(self._append_group_name("Vmin"),
                                                             self._append_group_name("Vmax"),
                                                             self._append_group_name("k0")))

        # self.beginIfBlock("{0} <= 0.0".format(self._append_group_name("k0doubleprime")))
        # self.beginIfBlock("{0} > 1.0".format(self._append_group_name("k0doubleprime")))
        # "k0doubleprime_window = (-k0doubleprime) * (1 - k0doubleprime)"
        self.addComputeGlobal(self._append_group_name("k0doubleprime_window"),
                              "(-{0}) * (1 - {0})".format(self._append_group_name("k0doubleprime")))

        self.beginIfBlock(self._append_group_name("k0doubleprime_window") + " >= 0.0")
        self.__calculate_simple_threshold_energy_and_effective_harmonic_constant()
        self.endBlock()

    def _lower_bound_calculate_threshold_energy_and_effective_harmonic_constant(self):
        self.__calculate_simple_threshold_energy_and_effective_harmonic_constant()
