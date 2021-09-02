from abc import ABC

from gamd.langevin.base_integrator import GroupBoostIntegrator
from simtk import unit as unit
from ..stage_integrator import BoostType
from ..stage_integrator import BoostMethod
from ..stage_integrator import ComputeType


class DualBoostIntegrator(GroupBoostIntegrator, ABC):
    def __init__(self, group, dt, ntcmdprep, ntcmd, ntebprep,
                 nteb, nstlim, ntave, sigma0p, sigma0d, collision_rate,
                 temperature, restart_filename):
        """
        Parameters
        ----------
        :param group:     The system group provided used by OpenMM for the Dihedral Energy and Forces.
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param nstlim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        :param sigma0p:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting. Applies to the total boost portion.
        :param sigma0d:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting. Applies to the dihedral boost portion.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """
        group_dict = {group: "Dihedral"}

        super(DualBoostIntegrator, self).__init__(group_dict,
                                                  BoostType.DUAL_TOTAL_DIHEDRAL,
                                                  BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL,
                                                  dt, ntcmdprep, ntcmd,
                                                  ntebprep, nteb, nstlim, ntave,
                                                  sigma0p, sigma0d,
                                                  collision_rate, temperature,
                                                  restart_filename)


class LowerBoundIntegrator(DualBoostIntegrator):
    def __init__(self, group, dt=2.0 * unit.femtoseconds, ntcmdprep=200000,
                 ntcmd=1000000, ntebprep=200000, nteb=1000000,
                 nstlim=3000000, ntave=50000,
                 sigma0p=6.0 * unit.kilocalories_per_mole,
                 sigma0d=6.0 * unit.kilocalories_per_mole,
                 collision_rate=1.0 / unit.picoseconds,
                 temperature=298.15 * unit.kelvin, restart_filename=None):
        """
        Parameters
        ----------
        :param group:     The system group provided used by OpenMM for the Dihedral Energy and Forces.
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param nstlim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        :param sigma0p:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting. Applies to the total boost portion.
        :param sigma0d:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting. Applies to the dihedral boost portion.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """
        self.__group = group
        super(LowerBoundIntegrator, self).__init__(group, dt, ntcmdprep, ntcmd, ntebprep,
                                                   nteb, nstlim, ntave, sigma0p, sigma0d,
                                                   collision_rate, temperature, restart_filename)
        
        # for i in range(self.getNumComputations()):
        #    print(self.getComputationStep(i))
        # exit()

    def _calculate_threshold_energy_and_effective_harmonic_constant(
            self, compute_type):
        super()._lower_bound_calculate_threshold_energy_and_effective_harmonic_constant(
            compute_type)


class UpperBoundIntegrator(DualBoostIntegrator):
    def __init__(self, group, dt=2.0 * unit.femtoseconds, ntcmdprep=200000,
                 ntcmd=1000000, ntebprep=200000, nteb=1000000,
                 nstlim=3000000, ntave=50000,
                 sigma0p=6.0 * unit.kilocalories_per_mole,
                 sigma0d=6.0 * unit.kilocalories_per_mole,
                 collision_rate=1.0 / unit.picoseconds,
                 temperature=298.15 * unit.kelvin, restart_filename=None):
        """
        Parameters
        ----------
        :param group:     The system group provided used by OpenMM for the Dihedral Energy and Forces.
        :param dt:        The Amount of time between each time step.
        :param ntcmdprep: The number of conventional MD steps for system equilibration.
        :param ntcmd:     The total number of conventional MD steps (including ntcmdprep). (must be a multiple of ntave)
        :param ntebprep:  The number of GaMD pre-equilibration steps.
        :param nteb:      The number of GaMD equilibration steps (including ntebprep). (must be a multiple of ntave)
        :param nstlim:    The total number of simulation steps.
        :param ntave:     The number of steps used to smooth the average and sigma of potential energy (corresponds to a
                          running average window size).
        :param sigma0p:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting. Applies to the total boost portion.
        :param sigma0d:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting. Applies to the dihedral boost portion.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """
        self.__group = group
        super(UpperBoundIntegrator, self).__init__(group, dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                                   nstlim, ntave, sigma0p, sigma0d,
                                                   collision_rate, temperature, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(
            self, compute_type):
        super()._upper_bound_calculate_threshold_energy_and_effective_harmonic_constant(
            compute_type)
