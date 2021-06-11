from abc import ABC

from gamd.langevin.base_integrator import GroupBoostIntegrator
from simtk import unit as unit
from ..stage_integrator import BoostType


class DihedralBoostIntegrator(GroupBoostIntegrator, ABC):
    def __init__(self, group, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, sigma0, collision_rate,
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
        :param sigma0:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """
        #self.__group = group
        #group_name = BoostType.DIHEDRAL
        group_dict = {group: "Dihedral"}
        total_boost = False

        super(DihedralBoostIntegrator, self).__init__(group_dict, total_boost, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim,
                                                      ntave, sigma0, collision_rate, temperature, restart_filename)


class LowerBoundIntegrator(DihedralBoostIntegrator):
    def __init__(self, group, dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000, ntebprep=200000, nteb=1000000,
                 nstlim=3000000, ntave=50000, sigma0=6.0 * unit.kilocalories_per_mole,
                 collision_rate=1.0 / unit.picoseconds, temperature=298.15 * unit.kelvin, restart_filename=None):
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
        :param sigma0:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """
        self.__group = 1
        super(LowerBoundIntegrator, self).__init__(group, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, sigma0,
                                                   collision_rate, temperature, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        super()._lower_bound_calculate_threshold_energy_and_effective_harmonic_constant()


class UpperBoundIntegrator(DihedralBoostIntegrator):
    def __init__(self, group, dt=2.0 * unit.femtoseconds, ntcmdprep=200000, ntcmd=1000000, ntebprep=200000, nteb=1000000,
                 nstlim=3000000, ntave=50000, sigma0=6.0 * unit.kilocalories_per_mole,
                 collision_rate=1.0 / unit.picoseconds, temperature=298.15 * unit.kelvin, restart_filename=None):
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
        :param sigma0:    The upper limit of the standard deviation of the potential boost that allows for
                          accurate reweighting.
        :param collision_rate:      Collision rate (gamma) compatible with 1/picoseconds, default: 1.0/unit.picoseconds
        :param temperature:         "Bath" temperature value compatible with units.kelvin, default: 298.15*unit.kelvin
        :param restart_filename:    The file name of the restart file.  (default=None indicates new simulation.)
        """
        self.__group = 1
        super(UpperBoundIntegrator, self).__init__(group, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, sigma0,
                                                      collision_rate, temperature, restart_filename)

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        super()._upper_bound_calculate_threshold_energy_and_effective_harmonic_constant()
