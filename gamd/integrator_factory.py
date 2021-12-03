"""
integrator_factory.py: Implements the GaMD integration method.

Portions copyright (c) 2021 University of Kansas
Authors: Matthew Copeland, Yinglong Miao
Contributors: Lane Votapka

"""

from simtk.unit import *

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


def set_all_forces_to_group(system):
    group = 1
    for force in system.getForces():
        force.setForceGroup(group)
    return group


def set_dihedral_group(system):
    return set_single_group(2, 'PeriodicTorsionForce', system)


def set_non_bonded_group(system):
    return set_single_group(1, 'NonbondedForce', system)


def set_single_group(group, name, system):
    for force in system.getForces():
        if force.__class__.__name__ == name:
            force.setForceGroup(group)
            break
    return group


def create_gamd_cmd_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    """
        This integrator is meant for use in generating a conventional MD baseline to compare against
        for the other integrators.

    :param system:
    :param temperature:
    :return:
    """
    group = set_dihedral_group(system)
    integrator = DihedralBoostLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                   ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                   ntave=ntave, temperature=temperature,
                                                   sigma0=0.0 * kilocalories_per_mole)
    result = ["", group, integrator]
    return result


def create_lower_total_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                                        sigma0=6.0 * kilocalories_per_mole):
    # The group is set, so that we can output the dihedral energy.  It doesn't impact calculations for total boost,
    # since we are utilizing the OpenMM provided variables with them not split out for total boost calculations.
    group = set_dihedral_group(system)
    integrator = TotalBoostLowerBoundIntegrator(dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                ntave=ntave, sigma0=sigma0, temperature=temperature)
    result = ["", group, integrator]
    return result


def create_upper_total_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                                        sigma0=6.0 * kilocalories_per_mole):
    # The group is set, so that we can output the dihedral energy.  It doesn't impact calculations for total boost,
    # since we are utilizing the OpenMM provided variables with them not split out for total boost calculations.
    group = set_dihedral_group(system)
    integrator = TotalBoostUpperBoundIntegrator(dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                ntave=ntave, sigma0=sigma0,
                                                temperature=temperature)
    result = ["", group, integrator]
    return result


def create_lower_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                                           sigma0=6.0 * kilocalories_per_mole):
    group = set_all_forces_to_group(system)
    integrator = DihedralBoostLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                                   nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0=sigma0,
                                                   temperature=temperature)
    result = ["", group, integrator]
    return result


def create_upper_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                                           sigma0=6.0 * kilocalories_per_mole):
    group = set_dihedral_group(system)
    integrator = DihedralBoostUpperBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                                   nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0=sigma0,
                                                   temperature=temperature)
    result = ["", group, integrator]
    return result


def create_lower_dual_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                                       sigma0p=6.0 * kilocalories_per_mole, sigma0d=6.0 * kilocalories_per_mole):
    group = set_dihedral_group(system)
    integrator = DualBoostLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                               nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0p=sigma0p,
                                               sigma0d=sigma0d, temperature=temperature)
    result = ["", group, integrator]
    return result


def create_upper_dual_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                                       sigma0p=6.0 * kilocalories_per_mole, sigma0d=6.0 * kilocalories_per_mole):
    group = set_dihedral_group(system)
    integrator = DualBoostUpperBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                               nteb=nteb, nstlim=nstlim, ntave=ntave,
                                               sigma0d=sigma0d,
                                               sigma0p=sigma0p, temperature=temperature)
    result = ["", group, integrator]
    return result


def create_lower_non_bonded_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                                             sigma0=6.0 * kilocalories_per_mole):
    group = set_non_bonded_group(system)
    integrator = NonBondedLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                               nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0=sigma0,
                                               temperature=temperature)
    result = ["", group, integrator]
    return result


def create_upper_non_bonded_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                                             sigma0=6.0 * kilocalories_per_mole):
    group = set_non_bonded_group(system)
    integrator = NonBondedUpperBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                               nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0=sigma0,
                                               temperature=temperature)
    result = ["", group, integrator]
    return result


def create_lower_dual_non_bonded_dihederal_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                            ntebprep, nteb, nstlim, ntave,
                                                            sigma0p=6.0 * kilocalories_per_mole,
                                                            sigma0d=6.0 * kilocalories_per_mole):
    nonbonded_group = set_non_bonded_group(system)
    dihedral_group = set_dihedral_group(system)
    integrator = DualNonBondedDihedralLowerIntegrator(nonbonded_group, dihedral_group, dt=dt, ntcmdprep=ntcmdprep,
                                                      ntcmd=ntcmd, ntebprep=ntebprep,
                                                      nteb=nteb, nstlim=nstlim, ntave=ntave, sigma0p=sigma0p,
                                                      sigma0d=sigma0d,
                                                      temperature=temperature)
    result = [nonbonded_group, dihedral_group, integrator]
    return result


def create_upper_dual_non_bonded_dihederal_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                            ntebprep, nteb, nstlim, ntave,
                                                            sigma0p=6.0 * kilocalories_per_mole,
                                                            sigma0d=6.0 * kilocalories_per_mole):
    nonbonded_group = set_non_bonded_group(system)
    dihedral_group = set_dihedral_group(system)
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
    def get_integrator(boost_type_str, system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave,
                       sigma0p=6.0 * kilocalories_per_mole, sigma0d=6.0 * kilocalories_per_mole):
        result = []
        first_boost_type = BoostType.TOTAL
        second_boost_type = BoostType.DIHEDRAL
        if boost_type_str == "gamd-cmd-base":
            result = create_gamd_cmd_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim,
                                                ntave)
        elif boost_type_str == "lower-total":
            result = create_lower_total_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                         ntebprep, nteb, nstlim, ntave, sigma0p)
        elif boost_type_str == "upper-total":
            result = create_upper_total_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                         ntebprep, nteb, nstlim, ntave, sigma0p)
        elif boost_type_str == "lower-dihedral":
            result = create_lower_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                            ntebprep, nteb, nstlim, ntave, sigma0p)
        elif boost_type_str == "upper-dihedral":
            result = create_upper_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                            ntebprep, nteb, nstlim, ntave, sigma0p)
        elif boost_type_str == "lower-dual":
            result = create_lower_dual_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                        ntebprep, nteb, nstlim, ntave, sigma0p, sigma0d)
        elif boost_type_str == "upper-dual":
            result = create_upper_dual_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                        ntebprep, nteb, nstlim, ntave, sigma0p, sigma0d)
        elif boost_type_str == "lower-nonbonded":
            result = create_lower_non_bonded_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                              ntebprep, nteb, nstlim, ntave, sigma0p)
            second_boost_type = BoostType.NON_BONDED
        elif boost_type_str == "upper-nonbonded":
            result = create_upper_non_bonded_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                              ntebprep, nteb, nstlim, ntave, sigma0p)
            second_boost_type = BoostType.NON_BONDED
        elif boost_type_str == "lower-dual-nonbonded-dihedral":
            result = create_upper_dual_non_bonded_dihederal_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                                             ntebprep, nteb, nstlim, ntave, sigma0p,
                                                                             sigma0d)
            first_boost_type = BoostType.NON_BONDED
            second_boost_type = BoostType.DIHEDRAL
        elif boost_type_str == "upper-dual-nonbonded-dihedral":
            result = create_upper_dual_non_bonded_dihederal_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                                             ntebprep, nteb, nstlim, ntave, sigma0p,
                                                                             sigma0d)
            first_boost_type = BoostType.NON_BONDED
            second_boost_type = BoostType.DIHEDRAL
        else:
            raise ValueError("Invalid boost_type_str passed to GamdIntegratorFactory.getIntegrator.")

        result.append(first_boost_type)
        result.append(second_boost_type)
        return result
