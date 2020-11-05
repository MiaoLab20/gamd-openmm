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
                                 "threshold_energy": -1E99,
                                 "boostPotential": 0.0,
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
        self.addUpdateContextState()
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
