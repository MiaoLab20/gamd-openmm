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
from gamd.langevin.base_integrator import GamdLangevinIntegrator


class TotalPotentialBoostIntegrator(GamdLangevinIntegrator, ABC):
    """ This class is an OpenMM Integrator for doing the total potential boost for
        Gaussian accelerated molecular dynamics.
    """

    def __init__(self,
                 dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000,
                 ntebprep=200000, nteb=1000000, ntslim=3000000, ntave=50000,
                 sigma0=6.0 * unit.kilocalories_per_mole,
                 collision_rate=1.0 / unit.picoseconds,
                 temperature=298.15 * unit.kelvin,
                 restart_filename=None):
        """
        Parameters
        ----------
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param ntslim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        :param sigma0:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """

        self.total_boost_global_variables = {"Vmax": -1E99, "Vmin": 1E99, "totalForceScalingFactor": 0, "Vavg": 0,
                                             "oldVavg": 0, "sigmaV": 0, "M2": 0,
                                             "wVavg": 0, "wVariance": 0, "k0": 0, "k0prime": 0, "k0doubleprime": 0,
                                             "currentPotentialEnergy": 0, "boosted_energy": 0, "sigma0": sigma0}
        self.total_boost_per_dof_variables = {"newx": 0, "coordinates": 0}
        self.debug_per_dof_variables = []
        # self.debug_per_dof_variables = ["x", "v", "f", "m"]
        self.debug_global_variables = ["dt", "energy", "energy0", "energy1", "energy2", "energy3", "energy4"]
        self.sigma0 = sigma0
        self.debuggingIsEnabled = True

        super(TotalPotentialBoostIntegrator, self).__init__(dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                                            ntslim, ntave, collision_rate, temperature,
                                                            restart_filename)
        # This makes sure that we get the value for the currentPotentialEnergy at the end of each simulation step.
        self.addComputeGlobal("currentPotentialEnergy", "energy")

        self.addComputePerDof("coordinates", "x")

    def get_starting_energy(self):
        return self.getGlobalVariableByName("starting_energy")

    def get_current_state(self):
        results = {"step": self.getGlobalVariableByName("stepCount")}

        for key, value in self.total_boost_global_variables.items():
            results[key] = self.getGlobalVariableByName(key)
        results["threshold_energy"] = self.getGlobalVariableByName("threshold_energy")
        results["boostPotential"] = self.getGlobalVariableByName("boostPotential")
        results["starting_energy"] = self.getGlobalVariableByName("starting_energy")

        
        #  results = {{ key:self.getGlobalVariableByName(key)} for key, value in self.total_boost_global_variables.items()}
        # garbage = {self.addPerDofVariable(key, value) for key, value in self.total_boost_per_dof_variables.items()}
        return results

    def _add_common_variables(self):
        garbage = {self.addGlobalVariable(key, value) for key, value in self.total_boost_global_variables.items()}
        garbage = {self.addPerDofVariable(key, value) for key, value in self.total_boost_per_dof_variables.items()}

        super(TotalPotentialBoostIntegrator, self)._add_common_variables()

    def _add_instructions_to_calculate_primary_boost_statistics(self):
        self.addComputeGlobal("Vmax", "max(energy, Vmax)")
        self.addComputeGlobal("Vmin", "min(energy, Vmin)")

    def _add_instructions_to_calculate_secondary_boost_statistics(self):
        #
        # The following calculations are used to keep a running average, rather than calculating
        # the average at the ntave % 0 step.
        #
        self.addComputeGlobal("oldVavg", "wVavg")
        self.addComputeGlobal("wVavg", "wVavg + (energy-wVavg)/windowCount")
        self.addComputeGlobal("M2", "M2 + (energy - oldVavg)*(energy - wVavg)")
        self.addComputeGlobal("wVariance", "select(windowCount-1,M2/(windowCount - 1),0)")

    def _update_potential_state_values_with_window_potential_state_values(self):
        # Update window variables
        self.addComputeGlobal("Vavg", "wVavg")
        self.addComputeGlobal("sigmaV", "sqrt(wVariance)")

        # reset variables
        self.addComputeGlobal("M2", "0")
        self.addComputeGlobal("wVavg", "0")
        self.addComputeGlobal("oldVavg", "0")
        self.addComputeGlobal("wVariance", "0")

    def _add_conventional_md_pre_calc_step(self):
        self.addComputePerDof("sigma", "sqrt(thermal_energy/m)")

    def _add_conventional_md_position_update_step(self):
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()

    def _add_conventional_md_velocity_update_step(self):
        #        self.addComputePerDof("v",
        #                              "v + 0.5*dt*fprime/m;fprime=f*((1.0-modify) + modify*(alpha/(alpha+threshold_energy-energy))^2); modify=step(threshold_energy-energy)")
        #
        #  The step function returns 0 for negative arguments and 1 for positive arguments.  Since the value of the
        #  threshold_energy defaults to -1E99 and thereshold_energy is not modified until the end of stage 2, the value
        #  of modify will always be zero in the equation above during the conventional md steps.
        #
        #  Given the above statements, we can reduce the above equation down the following form:
        #
        self.addComputePerDof("v", "v + 0.5*dt*f/m")
        self.addConstrainVelocities()

    def _add_conventional_md_stochastic_velocity_update_step(self):
        self.addComputePerDof("v", "current_velocity_component * v + random_velocity_component * sigma * gaussian")
        self.addConstrainVelocities()

    def _add_gamd_pre_calc_step(self):
        self.addComputePerDof("sigma", "sqrt(thermal_energy/m)")
        #
        # We don't actually apply the boost potential to the energy value, since energy is a read only value.
        #
        self.addComputeGlobal("boostPotential", "0.5 * k0 * (threshold_energy - energy)^2/(Vmax-Vmin)")
        self.addComputeGlobal("boosted_energy", "energy + boostPotential")

        self.beginIfBlock("boosted_energy >= threshold_energy")

        self.addComputeGlobal("boostPotential", "0")
        self.addComputeGlobal("boosted_energy", "energy")

        self.endBlock()

    def _add_gamd_boost_calculations_step(self):
        self.addComputeGlobal("totalForceScalingFactor", "1.0-((k0 * (threshold_energy - energy))/(Vmax - Vmin))")
        self.beginIfBlock("boosted_energy >= threshold_energy")
        self.addComputeGlobal("totalForceScalingFactor", "1.0")
        self.endBlock()

    def _add_gamd_position_update_step(self):
        self.addComputePerDof("x", "x + 0.5*dt*v")
        self.addComputePerDof("newx", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v + 2*(x-newx)/dt")
        self.addConstrainVelocities()

    def _add_gamd_velocity_update_step(self):
        self.addComputePerDof("v", "v + 0.5*(dt*f*totalForceScalingFactor/m)")
        self.addConstrainVelocities()

    def _add_gamd_stochastic_velocity_update_step(self):
        self.addComputePerDof("v", "current_velocity_component*v + random_velocity_component*sigma*gaussian")
        self.addConstrainVelocities()

    def get_current_potential_energy(self):
        return self.getGlobalVariableByName("currentPotentialEnergy")

    def get_force_scaling_factor(self):
        return self.getGlobalVariableByName("totalForceScalingFactor")

    def get_boost_potential(self):
        return self.getGlobalVariableByName("boostPotential")

    def get_boosted_total_energy(self):
        return self.getGlobalVariableByName("boosted_energy")

    def get_dihedral_boost(self):
        return 0.0

    #
    # Debugging Methods
    #

    def _add_debug(self):

        """
            This method will save off a copy of all of the variables in use.  Obviously, you don't
            really want to leave these statements in for a normal production run, since they will generate
            a lot of data that will slow down the simulation.

            It's recommended to put this call in an IfBlock of the particular step time frame you are interested
            in performing the debugging, rather than for an entire simulation.
        :return: none
        """

        if self.debuggingIsEnabled:
            garbage = {self._save_global_debug(key) for key, value in self.total_boost_global_variables.items()}
            garbage = [self._save_global_debug(name) for name in self.debug_global_variables]
            # garbage = {self._save_per_dof_debug(key) for key, value in self.total_boost_per_dof_variables.items()}
            garbage = [self._save_per_dof_debug(name) for name in self.debug_per_dof_variables]

            super(TotalPotentialBoostIntegrator, self)._add_debug()

    def get_debug_step(self, counter):
        """
            This method will retrieve all of the debugging values for a particular debug counter
            value.

        :param counter: an integer
        :return: a dictionary containing the global and per dof debug values for the simulation at that point.
        """
        results = super(TotalPotentialBoostIntegrator, self).get_debug_step(counter)
        #results.update({k: self._get_global_debug_value(counter, k) for (k, v) in self.total_boost_global_variables.items()})

        results.update(self._get_debug_values_as_dictionary(self.total_boost_global_variables, counter,
                                                            self._get_global_debug_value))
        # results.update(self._get_debug_values_as_dictionary(self.total_boost_per_dof_variables, counter,
        #                                                    self._get_per_dof_debug_value))

        for name in self.debug_per_dof_variables:
            results[str(counter) + "_" + name] = self._get_per_dof_debug_value(counter, name)

        for name in self.debug_global_variables:
            results[str(counter) + "_" + name] = self._get_global_debug_value(counter, name)

        return results

#    def _get_global_debug_value(self, counter, name):
#        super(TotalPotentialBoostIntegrator, self)._get_global_debug_value(counter, name)

    def get_debugging_information(self):
        """
            This method will retrieve all debugging values for all debug counters.  NOTE:  This data structure
            can be quite large.
        :return: a dictionary containing the global and per dof debug values for all debug time points.
        """
        debug_information = {}
        for count in range(0, self.debug_counter):
            debug_information[count] = self.get_debug_step(count)

        return debug_information


# ============================================================================================
# Integrator Classes
# ============================================================================================

# ============================================================================================
# Total Boost Integrator Classes
# ============================================================================================


class LowerBoundIntegrator(TotalPotentialBoostIntegrator):

    def __init__(self, dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000, ntebprep=200000, nteb=1000000,
                 ntslim=3000000, ntave=50000, sigma0=6.0 * unit.kilocalories_per_mole,
                 collision_rate=1.0 / unit.picoseconds, temperature=298.15 * unit.kelvin, restart_filename=None):
        """
        Parameters
        ----------
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param ntslim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        :param sigma0:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """

        super(LowerBoundIntegrator, self).__init__(dt, ntcmdprep, ntcmd,
                                                   ntebprep, nteb, ntslim,
                                                   ntave, sigma0, collision_rate,
                                                   temperature, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        self.addComputeGlobal("threshold_energy", "Vmax")
        self.addComputeGlobal("k0prime", "(sigma0/sigmaV) * (Vmax - Vmin)/(Vmax - Vavg)")
        self.addComputeGlobal("k0", "min(1.0, k0prime); ")


class UpperBoundIntegrator(TotalPotentialBoostIntegrator):

    def __init__(self, dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000,
                 ntebprep=200000, nteb=1000000, ntslim=3000000,
                 ntave=50000,
                 collision_rate=1.0 / unit.picoseconds,
                 temperature=298.15 * unit.kelvin,
                 sigma0=6.0 * unit.kilocalories_per_mole, restart_filename=None):
        """
        Parameters
        ----------
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param ntslim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        :param sigma0:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """

        super(UpperBoundIntegrator, self).__init__(dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                                   ntslim, ntave, collision_rate,
                                                   temperature, sigma0, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        self.addComputeGlobal("k0", "1.0")
        self.addComputeGlobal("k0doubleprime", "1-(sigma0/sigmaV) * (Vmax - Vmin)/(Vavg - Vmin)")
        #
        # We don't have an else statement in OpenMM, so we are going to set the value assuming that k0doubleprime
        # is between 0 and 1.  Then, we will override those values, if it is between them.
        #
        self.addComputeGlobal("k0", "k0doubleprime")
        self.addComputeGlobal("threshold_energy", "Vmin + (Vmax - Vmin)/k0")

        self.beginIfBlock("k0doubleprime <= 0.0")
        self.beginIfBlock("k0doubleprime > 1.0")
        self.addComputeGlobal("threshold_energy", "Vmax")
        self.addComputeGlobal("k0prime", "sigma0")

        self.endBlock()
        self.endBlock()
