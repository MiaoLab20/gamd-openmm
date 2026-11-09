"""
test_boost_force_update.py

Check that the velocity update of every boost type applies the force that
GaMD prescribes, using a small system whose forces live in group 0 plus the
boosted groups. The expected forces are built from per-group forces queried
from the Context, independently of how the integrator expresses the update.
"""

import numpy as np
import openmm
import openmm.unit as unit
import pytest

from gamd import integrator_factory as factory
from gamd.stage_integrator import BoostMethod

STAGE_ARGS = dict(dt=0.5 * unit.femtoseconds, ntcmdprep=20, ntcmd=40,
                  ntebprep=20, nteb=40, nstlim=10 ** 6, ntave=20)
PRODUCTION_STEP = 82
MEASURED_STEPS = 15

KJ = unit.kilojoule_per_mole
NM = unit.nanometer


def make_chain_system(n=8):
    """A bead chain with bonds, angles, torsions and Lennard-Jones/Coulomb."""
    system = openmm.System()
    for _ in range(n):
        system.addParticle(12.0)
    bonds = openmm.HarmonicBondForce()
    angles = openmm.HarmonicAngleForce()
    torsions = openmm.PeriodicTorsionForce()
    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    for i in range(n):
        nonbonded.addParticle(0.3 * (-1) ** i, 0.3, 0.5)
    for i in range(n - 1):
        bonds.addBond(i, i + 1, 0.15, 2.0e5)
    for i in range(n - 2):
        angles.addAngle(i, i + 1, i + 2, 1.9, 300.0)
    for i in range(n - 3):
        torsions.addTorsion(i, i + 1, i + 2, i + 3, 3, 0.0, 6.0)
        torsions.addTorsion(i, i + 1, i + 2, i + 3, 1, 0.5, 4.0)
    for i, j in [(i, j) for i in range(n) for j in range(i + 1, n)
                 if j - i <= 2]:
        nonbonded.addException(i, j, 0.0, 0.3, 0.0)
    for force in (bonds, angles, torsions, nonbonded):
        system.addForce(force)
    positions = [openmm.Vec3(0.15 * i, 0.1 * (i % 2), 0.05 * ((i // 2) % 2))
                 for i in range(n)]
    return system, positions


def group_forces(context, group):
    state = context.getState(getForces=True, groups={group})
    return state.getForces(asNumpy=True).value_in_unit(KJ / NM)


CASES = {
    "total": (factory.create_lower_total_boost_integrator, {}),
    "dihedral": (factory.create_lower_dihedral_boost_integrator,
                 {2: "Dihedral"}),
    "dual_total_dihedral": (factory.create_lower_dual_boost_integrator,
                            {2: "Dihedral"}),
    "non_bonded": (factory.create_lower_non_bonded_boost_integrator,
                   {1: "NonBonded"}),
    "dual_non_bonded_dihedral": (
        factory.create_lower_dual_non_bonded_dihederal_boost_integrator,
        {1: "NonBonded", 2: "Dihedral"}),
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_boosted_velocity_update_matches_gamd_forces(case):
    create, boosted_groups = CASES[case]
    system, positions = make_chain_system()
    integrator = create(system, 300 * unit.kelvin, **STAGE_ARGS)[2]
    integrator.setRandomNumberSeed(1234)
    context = openmm.Context(system, integrator,
                             openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    context.setVelocitiesToTemperature(300 * unit.kelvin, 5)
    integrator.step(PRODUCTION_STEP)
    assert integrator.getGlobalVariableByName("stage") >= 4

    # Turn off the stochastic kick so the update is deterministic.
    integrator.setGlobalVariableByName("noisescale", 0.0)
    integrator.setGlobalVariableByName("thermal_energy", 0.0)

    masses = np.array([system.getParticleMass(i).value_in_unit(unit.dalton)
                       for i in range(system.getNumParticles())])[:, None]
    boost_seen = False
    for _ in range(MEASURED_STEPS):
        factors = check_one_step(context, integrator, masses, boosted_groups)
        boost_seen = boost_seen or any(abs(k - 1.0) > 1e-3 for k in factors.values())
    # The boost must actually be on, otherwise this test proves nothing.
    assert boost_seen


def check_one_step(context, integrator, masses, boosted_groups):
    v0 = context.getState(getVelocities=True).getVelocities(
        asNumpy=True).value_in_unit(NM / unit.picosecond)
    forces = {g: group_forces(context, g) for g in [0] + list(boosted_groups)}
    total_force = context.getState(getForces=True).getForces(
        asNumpy=True).value_in_unit(KJ / NM)
    vscale = integrator.getGlobalVariableByName("vscale")
    fscale = integrator.getGlobalVariableByName("fscale")

    integrator.step(1)

    v1 = context.getState(getVelocities=True).getVelocities(
        asNumpy=True).value_in_unit(NM / unit.picosecond)
    applied = (v1 - vscale * v0) * masses / fscale

    factors = {name: integrator.getGlobalVariableByName(
        "ForceScalingFactor_" + name) for name in boosted_groups.values()}
    method = integrator._boost_method
    if method == BoostMethod.TOTAL:
        k_total = integrator.getGlobalVariableByName("ForceScalingFactor_Total")
        expected = k_total * total_force
        factors["Total"] = k_total
    else:
        expected = forces[0].copy()
        for group, name in boosted_groups.items():
            expected += factors[name] * forces[group]
        if method == BoostMethod.DUAL_DEPENDENT_GROUP_TOTAL:
            k_total = integrator.getGlobalVariableByName(
                "ForceScalingFactor_Total")
            expected *= k_total
            factors["Total"] = k_total

    np.testing.assert_allclose(applied, expected, rtol=1e-6,
                               atol=1e-6 * np.abs(expected).max())
    return factors


def test_dual_boost_total_term_shares_group_with_energy():
    """
    The total-boost term must use "f" (all groups) rather than "f0" so that
    OpenMM can evaluate it together with the "energy" request in one pass.
    """
    system, _ = make_chain_system()
    integrator = factory.create_lower_dual_boost_integrator(
        system, 300 * unit.kelvin, **STAGE_ARGS)[2]
    per_dof_expressions = []
    for i in range(integrator.getNumComputations()):
        kind, variable, expression = integrator.getComputationStep(i)
        if kind == openmm.CustomIntegrator.ComputePerDof and variable == "v":
            per_dof_expressions.append(expression)
    assert not any("f0" in expr for expr in per_dof_expressions)
