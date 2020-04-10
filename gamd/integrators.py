"""
gamd.py: Implements the GaMD integration method.

Portions copyright (c) 2020 University of Kansas
Authors: Matthew Copeland, Yinglong Miao
Contributors: Lane Votapka

"""

from __future__ import absolute_import

__author__ = "Matthew Copeland"
__version__ = "1.0"

from simtk.openmm import CustomIntegrator
from simtk import unit as unit
from abc import ABCMeta, ABC
from abc import abstractmethod
import numpy


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
                 ntebprep=200000, nteb=1000000, ntslim=3000000, ntave=50000):
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
        self.stage_3_start = ntcmd + 1
        self.stage_3_end = ntcmd + ntebprep
        self.stage_4_start = ntcmd + ntebprep + 1
        self.stage_4_end = ntcmd + nteb
        self.stage_5_start = ntcmd + nteb + 1
        self.stage_5_end = ntslim

        self.dt = dt
        self.ntcmdprep = ntcmdprep
        self.ntcmd = ntcmd
        self.ntebprep = ntebprep
        self.nteb = nteb
        self.ntslim = ntslim
        self.ntave = ntave

        self.addGlobalVariable("stepCount", 0)
        self.addGlobalVariable("windowCount", 0)
        self.addGlobalVariable("stage", -1)
        self.addComputeGlobal("stepCount", "stepCount+1")

        self._add_common_variables()
        self._run_precalculations()
        self.addUpdateContextState()

        self._add_stage_one_instructions()
        self._add_stage_two_instructions()
        self._add_stage_three_instructions()
        self._add_stage_four_instructions()
        self._add_stage_five_instructions()

    def _add_stage_one_instructions(self):
        self.beginIfBlock("stepCount <= " + str(self.stage_1_end))
        # -------------------------------
        self.addComputeGlobal("stage", "1")
        self._add_conventional_md_instructions()
        # -------------------------------
        self.endBlock()

    def _add_stage_three_instructions(self):
        self.beginIfBlock("stepCount >= " + str(self.stage_3_start))
        self.beginIfBlock("stepCount <= " + str(self.stage_3_end))
        # -------------------------------
        self.addComputeGlobal("stage", "3")
        self._add_apply_boost_potential_calculations()
        # -------------------------------
        self.endBlock()
        self.endBlock()

    def _add_stage_five_instructions(self):
        self.beginIfBlock("stepCount >= " + str(self.stage_5_start))
        self.beginIfBlock("stepCount <= " + str(self.stage_5_end))
        # -------------------------------
        self.addComputeGlobal("stage", "5")
        self._add_apply_boost_potential_calculations()
        # -------------------------------
        self.endBlock()
        self.endBlock()

    def _add_stage_two_instructions(self):
        self.beginIfBlock("stepCount >= " + str(self.stage_2_start))
        self.beginIfBlock("stepCount <= " + str(self.stage_2_end))

        # -------------------------------
        self.addComputeGlobal("stage", "2")
        self.addComputeGlobal("windowCount", "windowCount + 1")

        self._add_conventional_md_instructions()

        # Check to see if we need to calculate boost parameters
        self.beginIfBlock("windowCount = " + str(self.ntave))
        self._add_boost_parameter_calculations()
        self.addComputeGlobal("windowCount", "0")
        self.endBlock()

        self.beginIfBlock("stepCount = " + str(self.stage_2_end))
        self._set_boost_parameters_for_future_steps()
        self.endBlock()

        # -------------------------------
        self.endBlock()
        self.endBlock()

    def _add_stage_four_instructions(self):
        self.beginIfBlock("stepCount = " + str(self.stage_4_start))
        self.addComputeGlobal("windowCount", "0")
        self.endBlock()

        self.beginIfBlock("stepCount >= " + str(self.stage_4_start))
        self.beginIfBlock("stepCount <= " + str(self.stage_4_end))
        # -------------------------------
        self.addComputeGlobal("stage", "4")
        self.addComputeGlobal("windowCount", "windowCount + 1")


        # Check to see if we need to calculate boost parameters
        self.beginIfBlock("windowCount = " + str(self.ntave))
        self._add_boost_parameter_calculations()
        self.addComputeGlobal("windowCount", "0")
        self.endBlock()

        self._set_boost_parameters_for_future_steps()
        self._add_apply_boost_potential_calculations()

        # -------------------------------
        self.endBlock()
        self.endBlock()

    def get_stage(self):
        return self.getGlobalVariableByName("stage")

    def get_step_count(self):
        return self.getGlobalVariableByName("stepCount")

    def get_window_count(self):
        return self.getGlobalVariableByName("windowCount")

    def get_total_simulation_steps(self):
        return self.ntslim

    def get_coordinates(self):
        return self.getPerDofVariableByName("coordinates")

    def create_positions_file(self, filename):
        positions = self.get_coordinates()
        with open(filename, 'w') as file:
            file.write("particle, x, y, z\n")
            for i in range(len(positions)):
                line = str(i) + ", " + str(positions[i][0]) + ", " + str(positions[i][1]) + ", " + str(positions[i][2])
                file.write(line + "\n")

    @abstractmethod
    def _add_common_variables(self):
        raise NotImplementedError("must implement _add_common_variables")

    @abstractmethod
    def _run_precalculations(self):
        raise NotImplementedError("must implement _run_precalcuations")

    @abstractmethod
    def _add_conventional_md_instructions(self):
        raise NotImplementedError("must implement _add_conventional_md_instructions")

    @abstractmethod
    def _add_boost_parameter_calculations(self):
        raise NotImplementedError("must implement _add_boost_parameter_calculations")

    @abstractmethod
    def _set_boost_parameters_for_future_steps(self):
        raise NotImplementedError("must implement _set_boost_parameters_for_future_steps")

    @abstractmethod
    def _add_apply_boost_potential_calculations(self):
        raise NotImplementedError("must implement _add_apply_boost_potential_calculations")

    @abstractmethod
    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        raise NotImplementedError("must implement _calculate_threshold_energy_and_effective_harmonic_constant")


class GamdLangevinIntegrator(GamdStageIntegrator, ABC):

    def __init__(self,
                 dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000,
                 ntebprep=200000, nteb=1000000, ntslim=3000000, ntave=50000,
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
         :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
         :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
         :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
         """

        self.collision_rate = collision_rate
        self.temperature = temperature
        self.restart_filename = restart_filename
        self.kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
        self.thermal_energy = self.kB * self.temperature  # kT
        self.current_velocity_component = numpy.exp(-self.collision_rate * dt)  # a
        self.random_velocity_component = numpy.sqrt(1 - numpy.exp(- 2 * self.collision_rate * dt))  # b

        #
        # We need to run our super classes constructor last, since it's going to execute our other methods, which
        # have dependencies on our variables above being setup.
        #
        super(GamdLangevinIntegrator, self).__init__(dt, ntcmdprep, ntcmd, ntebprep, nteb, ntslim, ntave)

    def _add_common_variables(self):
        #
        #  totalEnergy = energy
        #
        self.addGlobalVariable("thermal_energy", self.thermal_energy)
        self.addGlobalVariable("current_velocity_component", self.current_velocity_component)
        self.addGlobalVariable("random_velocity_component", self.random_velocity_component)
        self.addGlobalVariable("collision_rate", self.collision_rate)
        self.addPerDofVariable("sigma", 0)

        self.addGlobalVariable("Vmax", -1E99)
        self.addGlobalVariable("Vmin", 1E99)


        # Unverified
        self.addGlobalVariable("Vavg", 0)
        self.addGlobalVariable("oldVavg", 0)
        self.addGlobalVariable("sigmaV", 0)

        self.addGlobalVariable("M2", 0)
        self.addGlobalVariable("wVmax", -1E99)
        self.addGlobalVariable("wVmin", 1E99)
        self.addGlobalVariable("wVavg", 0)
        self.addGlobalVariable("wVariance", 0)




class GamdTotalPotentialBoostLangevinIntegrator(GamdLangevinIntegrator, ABC):
    """ This class is an OpenMM Integrator for doing the total potential boost for
        Gaussian accelerated molecular dynamics.
    """

    def __init__(self,
                 dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000,
                 ntebprep=200000, nteb=1000000,ntslim=3000000, ntave=50000,
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

        self.sigma0 = sigma0
        super(GamdTotalPotentialBoostLangevinIntegrator, self).__init__(dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                                                        ntslim, ntave, collision_rate, temperature,
                                                                        restart_filename)

    def _add_common_variables(self):
        super(GamdTotalPotentialBoostLangevinIntegrator, self)._add_common_variables()

    def _run_precalculations(self):
        self.addComputeGlobal("Vmax", "max(Vmax, energy)")
        self.addComputeGlobal("Vmin", "min(Vmin, energy)")

    def _add_conventional_md_instructions(self):
        pass

    def _add_boost_parameter_calculations(self):
        pass

    def _set_boost_parameters_for_future_steps(self):
        pass

    def _add_apply_boost_potential_calculations(self):
        pass


class GamdDihedralBoostLangevinIntegrator(GamdLangevinIntegrator, ABC):
    """ This class is an OpenMM Integrator for doing the Dihedral Boost for
        Gaussian accelerated molecular dynamics.
    """

    def __init__(self, group,
                 dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000, ntebprep=200000, nteb=1000000,
                 ntslim=3000000, ntave=50000, sigma0=6.0 * unit.kilocalories_per_mole,
                 collision_rate=1.0 / unit.picoseconds, temperature=298.15 * unit.kelvin, restart_filename=None):
        """
        Parameters
        ----------
        :param group      The set of particles to boost.
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

        self.group = group
        self.sigma0 = sigma0
        super(GamdDihedralBoostLangevinIntegrator, self).__init__(dt, ntcmdprep, ntcmd, ntebprep, nteb, ntslim, ntave,
                                                                  collision_rate, temperature, restart_filename)

    def _add_common_variables(self):
        super(GamdDihedralBoostLangevinIntegrator, self)._add_common_variables()

    def _add_conventional_md_instructions(self):
        pass

    def _add_boost_parameter_calculations(self):
        pass

    def _set_boost_parameters_for_future_steps(self):
        pass

    def _add_apply_boost_potential_calculations(self):
        pass


class GamdDualBoostLangevinIntegrator(GamdLangevinIntegrator, ABC):
    """ This class is an OpenMM Integrator for doing the boost on both the total potential energy and the
        dihedral energy for Gaussian accelerated molecular dynamics.
    """

    def __init__(self, group, dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000, ntebprep=200000,
                 nteb=1000000, ntslim=3000000, ntave=50000, sigma0D=6.0 * unit.kilocalories_per_mole,
                 sigma0P=6.0 * unit.kilocalories_per_mole, collision_rate=1.0 / unit.picoseconds,
                 temperature=298.15 * unit.kelvin, restart_filename=None):
        """
        Parameters
        ----------
        :param group      The set of particles to boost.
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param ntslim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        :param sigma0D:   The upper limit of the standard deviation of the dihedral potential boost that allows for
                          accurate reweighting.
        :param sigma0P:   The upper limit of the standard deviation of the total potential boost that allows for
                          accurate reweighting.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """

        self.group = group
        self.sigma0D = sigma0D
        self.sigma0P = sigma0P
        super(GamdDualBoostLangevinIntegrator, self).__init__(dt, ntcmdprep, ntcmd, ntebprep, nteb, ntslim, ntave,
                                                              collision_rate, temperature, restart_filename)

    def _add_common_variables(self):
        super(GamdDualBoostLangevinIntegrator, self)._add_common_variables()

    def _add_conventional_md_instructions(self):
        pass

    def _add_boost_parameter_calculations(self):
        pass

    def _set_boost_parameters_for_future_steps(self):
        pass

    def _add_apply_boost_potential_calculations(self):
        pass


# ============================================================================================
# Integrator Classes
# ============================================================================================

# ============================================================================================
# Total Boost Integrator Classes
# ============================================================================================

class GamdTotalPotentialBoostLangevinLowerBoundIntegrator(GamdTotalPotentialBoostLangevinIntegrator):

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

        super(GamdTotalPotentialBoostLangevinLowerBoundIntegrator, self).__init__(dt, ntcmdprep, ntcmd,
                                                                                  ntebprep, nteb, ntslim,
                                                                                  ntave,  sigma0, collision_rate,
                                                                                  temperature, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        pass


class GamdTotalPotentialBoostLangevinUpperBoundIntegrator(GamdTotalPotentialBoostLangevinIntegrator):

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

        super(GamdTotalPotentialBoostLangevinUpperBoundIntegrator, self).__init__(dt, ntcmdprep, ntcmd,ntebprep, nteb,
                                                                                  ntslim, ntave, collision_rate,
                                                                                  temperature, sigma0, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        pass


# ============================================================================================
# Dihedral Boost Integrator Classes
# ============================================================================================


class GamdDihedralBoostLangevinLowerBoundIntegrator(GamdDihedralBoostLangevinIntegrator):

    def __init__(self,  group, dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000, ntebprep=200000,
                 nteb=1000000, ntslim=3000000, ntave=50000, sigma0=6.0 * unit.kilocalories_per_mole,
                 collision_rate=1.0 / unit.picoseconds, temperature=298.15 * unit.kelvin, restart_filename=None):
        """
        Parameters
        ----------
        :param group      The set of particles to boost.
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

        super(GamdDihedralBoostLangevinLowerBoundIntegrator, self).__init__(group, dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                                                            ntslim, ntave, sigma0, collision_rate,
                                                                            temperature, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        pass


class GamdDihedralBoostLangevinUpperBoundIntegrator(GamdDihedralBoostLangevinIntegrator):

    def __init__(self,  group, dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000, ntebprep=200000,
                 nteb=1000000, ntslim=3000000, ntave=50000, sigma0=6.0 * unit.kilocalories_per_mole,
                 collision_rate=1.0 / unit.picoseconds, temperature=298.15 * unit.kelvin, restart_filename=None):
        """
        Parameters
        ----------
        :param group      The set of particles to boost.
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

        super(GamdDihedralBoostLangevinUpperBoundIntegrator, self).__init__(group, dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                                                            ntslim, ntave, sigma0, collision_rate,
                                                                            temperature, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        pass


# ============================================================================================
# Dual Boost Integrator Classes
# ============================================================================================


class GamdDualBoostLangevinLowerBoundIntegrator(GamdDualBoostLangevinIntegrator):

    def __init__(self, group, dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000, ntebprep=200000,
                 nteb=1000000, ntslim=3000000, ntave=50000, sigma0D=6.0 * unit.kilocalories_per_mole,
                 sigma0P=6.0 * unit.kilocalories_per_mole, collision_rate=1.0 / unit.picoseconds,
                 temperature=298.15 * unit.kelvin, restart_filename=None):
        """
        Parameters
        ----------
        :param group      The set of particles to boost.
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param ntslim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        :param sigma0D:   The upper limit of the standard deviation of the dihedral potential boost that allows for
                          accurate reweighting.
        :param sigma0P:   The upper limit of the standard deviation of the total potential boost that allows for
                          accurate reweighting.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """
        super(GamdDualBoostLangevinLowerBoundIntegrator, self).__init__(group, dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                                                        ntslim, ntave, collision_rate, temperature,
                                                                        sigma0D, sigma0P, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        pass


class GamdDualBoostLangevinUpperBoundIntegrator(GamdDualBoostLangevinIntegrator):

    def __init__(self, group, dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000, ntebprep=200000,
                 nteb=1000000, ntslim=3000000, ntave=50000, sigma0D=6.0 * unit.kilocalories_per_mole,
                 sigma0P=6.0 * unit.kilocalories_per_mole, collision_rate=1.0 / unit.picoseconds,
                 temperature=298.15 * unit.kelvin, restart_filename=None):
        """
        Parameters
        ----------
        :param group      The set of particles to boost.
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param ntslim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        :param sigma0D:   The upper limit of the standard deviation of the dihedral potential boost that allows for
                          accurate reweighting.
        :param sigma0P:   The upper limit of the standard deviation of the total potential boost that allows for
                          accurate reweighting.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """
        super(GamdDualBoostLangevinUpperBoundIntegrator, self).__init__(group, dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                                                        ntslim, ntave, collision_rate, temperature,
                                                                        sigma0D, sigma0P, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        pass
