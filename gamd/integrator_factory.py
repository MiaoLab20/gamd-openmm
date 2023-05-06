"""
integrator_factory.py: Implements the GaMD integration method.

Portions copyright (c) 2021 University of Kansas
Authors: Matthew Copeland, Yinglong Miao
Contributors: Lane Votapka

"""

import openmm.unit as unit

from gamd.langevin.dihedral_boost_integrators import LowerBoundIntegrator as DihedralBoostLowerBoundIntegrator
from gamd.langevin.dihedral_boost_integrators import UpperBoundIntegrator as DihedralBoostUpperBoundIntegrator
from gamd.langevin.dual_boost_integrators import LowerBoundIntegrator as DualBoostLowerBoundIntegrator
from gamd.langevin.dual_boost_integrators import UpperBoundIntegrator as DualBoostUpperBoundIntegrator
from gamd.langevin.dual_non_bonded_dihedral_boost_integrators import \
    LowerBoundIntegrator as DualNonBondedDihedralLowerIntegrator
from gamd.langevin.dual_non_bonded_dihedral_boost_integrators import \
    UpperBoundIntegrator as DualNonBondedDihedralUpperIntegrator
from gamd.langevin.non_bonded_boost_integrators import LowerBoundIntegrator as NonBondedLowerBoundIntegrator
from gamd.langevin.non_bonded_boost_integrators import UpperBoundIntegrator as NonBondedUpperBoundIntegrator
from gamd.langevin.total_boost_integrators import LowerBoundIntegrator as TotalBoostLowerBoundIntegrator
from gamd.langevin.total_boost_integrators import UpperBoundIntegrator as TotalBoostUpperBoundIntegrator
from gamd.stage_integrator import BoostType


def print_force_group_information(system):
    for force in system.getForces():
        print("Force Name, Group:  ", force.__class__.__name__,
              ", ", force.getForceGroup())


def set_all_forces_to_group(system, group):
    for force in system.getForces():
        force.setForceGroup(group)


def set_dihedral_group(config, system):
    return set_single_group(2, ['PeriodicTorsionForce', 'CMAPTorsionForce'], system)


def set_non_bonded_group(config, system):
    return set_single_group(1, ['NonbondedForce', 'CustomNonbondedForce'], system)


def set_single_group(group, name_list, system):
    for force in system.getForces():
        if force.__class__.__name__ in name_list:
            # print(force.__class__.__name__)
            force.setForceGroup(group)
    return group

def create_gamd_cmd_integrator(config, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    """
        This integrator is meant for use in generating a conventional MD baseline to compare against
        for the other integrators.

    :param system:
    :param temperature:
    :return:
    """
    group = set_dihedral_group(config,system)
    integrator = DihedralBoostLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                   ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                   ntave=ntave, temperature=temperature,
                                                   sigma0=0.0 * unit.kilocalories_per_mole)
    result = ["", group, integrator]
    return result


def create_lower_total_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim,
                                        ntave,
                                        sigma0=6.0 * unit.kilocalories_per_mole):
    # The group is set, so that we can output the dihedral energy.  It doesn't impact calculations for total boost,
    # since we are utilizing the OpenMM provided variables with them not split out for total boost calculations.
    group = set_dihedral_group(config,system)
    integrator = TotalBoostLowerBoundIntegrator(dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                ntave=ntave, sigma0=sigma0, temperature=temperature)
    result = ["", group, integrator]
    return result


def create_upper_total_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim,
                                        ntave,
                                        sigma0=6.0 * unit.kilocalories_per_mole):
    # The group is set, so that we can output the dihedral energy.  It doesn't impact calculations for total boost,
    # since we are utilizing the OpenMM provided variables with them not split out for total boost calculations.
    group = set_dihedral_group(config,system)
    integrator = TotalBoostUpperBoundIntegrator(dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                ntave=ntave, sigma0=sigma0,
                                                temperature=temperature)
    result = ["", group, integrator]
    return result


def create_lower_dihedral_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim,
                                           ntave,
                                           sigma0=6.0 * unit.kilocalories_per_mole):
    group = set_dihedral_group(config,system)
    integrator = DihedralBoostLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                                   nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0=sigma0,
                                                   temperature=temperature)
    result = ["", group, integrator]
    return result


def create_upper_dihedral_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim,
                                           ntave,
                                           sigma0=6.0 * unit.kilocalories_per_mole):
    group = set_dihedral_group(config,system)
    integrator = DihedralBoostUpperBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                                   nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0=sigma0,
                                                   temperature=temperature)
    result = ["", group, integrator]
    return result


def create_lower_dual_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                                       sigma0p=6.0 * unit.kilocalories_per_mole,
                                       sigma0d=6.0 * unit.kilocalories_per_mole):
    group = set_dihedral_group(config,system)
    integrator = DualBoostLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                               nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0p=sigma0p,
                                               sigma0d=sigma0d, temperature=temperature)
    result = ["", group, integrator]
    return result


def create_upper_dual_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                                       sigma0p=6.0 * unit.kilocalories_per_mole,
                                       sigma0d=6.0 * unit.kilocalories_per_mole):
    group = set_dihedral_group(config,system)
    integrator = DualBoostUpperBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                               nteb=nteb, nstlim=nstlim, ntave=ntave,
                                               sigma0d=sigma0d,
                                               sigma0p=sigma0p, temperature=temperature)
    result = ["", group, integrator]
    return result


def create_lower_non_bonded_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim,
                                             ntave,
                                             sigma0=6.0 * unit.kilocalories_per_mole):
    group = set_non_bonded_group(config,system)
    integrator = NonBondedLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                               nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0=sigma0,
                                               temperature=temperature)
    result = ["", group, integrator]
    return result


def create_upper_non_bonded_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim,
                                             ntave,
                                             sigma0=6.0 * unit.kilocalories_per_mole):
    group = set_non_bonded_group(config,system)
    integrator = NonBondedUpperBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                               nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0=sigma0,
                                               temperature=temperature)
    result = ["", group, integrator]
    return result


def create_lower_dual_non_bonded_dihederal_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd,
                                                            ntebprep, nteb, nstlim, ntave,
                                                            sigma0p=6.0 * unit.kilocalories_per_mole,
                                                            sigma0d=6.0 * unit.kilocalories_per_mole):
    nonbonded_group = set_non_bonded_group(config,system)
    dihedral_group = set_dihedral_group(config,system)
    integrator = DualNonBondedDihedralLowerIntegrator(nonbonded_group, dihedral_group, dt=dt, ntcmdprep=ntcmdprep,
                                                      ntcmd=ntcmd, ntebprep=ntebprep,
                                                      nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0p=sigma0p,
                                                      sigma0d=sigma0d,
                                                      temperature=temperature)
    result = [nonbonded_group, dihedral_group, integrator]
    return result


def create_upper_dual_non_bonded_dihederal_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd,
                                                            ntebprep, nteb, nstlim, ntave,
                                                            sigma0p=6.0 * unit.kilocalories_per_mole,
                                                            sigma0d=6.0 * unit.kilocalories_per_mole):
    nonbonded_group = set_non_bonded_group(config,system)
    dihedral_group = set_dihedral_group(config,system)
    integrator = DualNonBondedDihedralUpperIntegrator(nonbonded_group, dihedral_group, dt=dt, ntcmdprep=ntcmdprep,
                                                      ntcmd=ntcmd, ntebprep=ntebprep,
                                                      nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0p=sigma0p,
                                                      sigma0d=sigma0d,
                                                      temperature=temperature)
    result = [nonbonded_group, dihedral_group, integrator]
    return result


class GamdIntegratorFactory:

    def __init__(self):
        pass

    @staticmethod
    def get_integrator(config, system):

        boost_type_str = config.integrator.boost_type
        temperature = config.temperature
        dt = config.integrator.dt
        ntcmdprep = config.integrator.number_of_steps.conventional_md_prep
        ntcmd = config.integrator.number_of_steps.conventional_md
        ntebprep = config.integrator.number_of_steps.gamd_equilibration_prep
        nteb = config.integrator.number_of_steps.gamd_equilibration
        nstlim = config.integrator.number_of_steps.total_simulation_length
        ntave = config.integrator.number_of_steps.averaging_window_interval
        sigma0p = config.integrator.sigma0.primary
        sigma0d = config.integrator.sigma0.secondary

        return GamdIntegratorFactory.get_integrator_by_parts(config, boost_type_str,
                                                             system, temperature, dt,
                                                             ntcmdprep, ntcmd, ntebprep, nteb,
                                                             nstlim, ntave, sigma0p, sigma0d)

    @staticmethod
    def get_integrator_by_parts(config, boost_type_str, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                nstlim, ntave,
                                sigma0p=6.0 * unit.kilocalories_per_mole, sigma0d=6.0 * unit.kilocalories_per_mole):
        set_all_forces_to_group(system,0)
        result = []
        first_boost_type = BoostType.TOTAL
        second_boost_type = BoostType.DIHEDRAL
        if boost_type_str == "gamd-cmd-base":
            result = create_gamd_cmd_integrator(config, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                                nstlim,
                                                ntave)
        elif boost_type_str == "lower-total":
            result = create_lower_total_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd,
                                                         ntebprep, nteb, nstlim, ntave, sigma0p)
        elif boost_type_str == "upper-total":
            result = create_upper_total_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd,
                                                         ntebprep, nteb, nstlim, ntave, sigma0p)
        elif boost_type_str == "lower-dihedral":
            result = create_lower_dihedral_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd,
                                                            ntebprep, nteb, nstlim, ntave, sigma0p)
        elif boost_type_str == "upper-dihedral":
            result = create_upper_dihedral_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd,
                                                            ntebprep, nteb, nstlim, ntave, sigma0p)
        elif boost_type_str == "lower-dual":
            result = create_lower_dual_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd,
                                                        ntebprep, nteb, nstlim, ntave, sigma0p, sigma0d)
        elif boost_type_str == "upper-dual":
            result = create_upper_dual_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd,
                                                        ntebprep, nteb, nstlim, ntave, sigma0p, sigma0d)
        elif boost_type_str == "lower-nonbonded":
            result = create_lower_non_bonded_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd,
                                                              ntebprep, nteb, nstlim, ntave, sigma0p)
            second_boost_type = BoostType.NON_BONDED
        elif boost_type_str == "upper-nonbonded":
            result = create_upper_non_bonded_boost_integrator(config, system, temperature, dt, ntcmdprep, ntcmd,
                                                              ntebprep, nteb, nstlim, ntave, sigma0p)
            second_boost_type = BoostType.NON_BONDED
        elif boost_type_str == "lower-dual-nonbonded-dihedral":
            result = create_lower_dual_non_bonded_dihederal_boost_integrator(config, system, temperature, dt, ntcmdprep,
                                                                             ntcmd,
                                                                             ntebprep, nteb, nstlim, ntave, sigma0p,
                                                                             sigma0d)
            first_boost_type = BoostType.NON_BONDED
            second_boost_type = BoostType.DIHEDRAL
        elif boost_type_str == "upper-dual-nonbonded-dihedral":
            result = create_upper_dual_non_bonded_dihederal_boost_integrator(config, system, temperature, dt, ntcmdprep,
                                                                             ntcmd,
                                                                             ntebprep, nteb, nstlim, ntave, sigma0p,
                                                                             sigma0d)
            first_boost_type = BoostType.NON_BONDED
            second_boost_type = BoostType.DIHEDRAL
        else:
            raise ValueError("Invalid boost_type_str passed to GamdIntegratorFactory.getIntegrator.")

        result.append(first_boost_type)
        result.append(second_boost_type)
        return result
