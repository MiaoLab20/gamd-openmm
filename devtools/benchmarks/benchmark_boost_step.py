"""
Time the production-stage step of each GaMD boost type.

Builds a solvated Trp-cage system (PDB 1L2Y, downloaded from RCSB if no PDB
is given) with PME, runs each integrator through its conventional MD and
equilibration stages, then times production steps. A plain Langevin
CustomIntegrator that requests forces and energy for all groups (one fused
force evaluation per step) is timed as a reference.

Usage:
    python benchmark_boost_step.py [--pdb file.pdb] [--platform CPU]
                                   [--steps 500] [--repeats 5] [--only dual]
"""

import argparse
import os
import statistics
import time
import urllib.request

import openmm
import openmm.app as app
import openmm.unit as unit

from gamd import integrator_factory as factory

STAGE_ARGS = dict(dt=2.0 * unit.femtoseconds, ntcmdprep=50, ntcmd=100,
                  ntebprep=50, nteb=100, nstlim=10 ** 9, ntave=50)
WARMUP_STEPS = 250  # past ntcmd + nteb, so we are in production

BOOSTS = {
    "total (lower)": factory.create_lower_total_boost_integrator,
    "dihedral (lower)": factory.create_lower_dihedral_boost_integrator,
    "dual total+dihedral (lower)": factory.create_lower_dual_boost_integrator,
    "non-bonded (lower)": factory.create_lower_non_bonded_boost_integrator,
    "dual non-bonded+dihedral (lower)":
        factory.create_lower_dual_non_bonded_dihederal_boost_integrator,
}


def build_system(pdb_path):
    if pdb_path is None:
        pdb_path = "1L2Y.pdb"
        if not os.path.exists(pdb_path):
            urllib.request.urlretrieve(
                "https://files.rcsb.org/download/1L2Y.pdb", pdb_path)
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    pdb = app.PDBFile(pdb_path)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, padding=1.0 * unit.nanometer)
    return forcefield, modeller


def reference_integrator(temperature):
    integrator = openmm.CustomIntegrator(STAGE_ARGS["dt"])
    integrator.addGlobalVariable("a", 0.0)
    integrator.addGlobalVariable("b", 0.0)
    integrator.addGlobalVariable("kT", (unit.MOLAR_GAS_CONSTANT_R * temperature))
    integrator.addGlobalVariable("e", 0.0)
    integrator.addPerDofVariable("x1", 0.0)
    integrator.addComputeGlobal("a", "exp(-dt)")
    integrator.addComputeGlobal("b", "sqrt(1-a*a)")
    integrator.addComputeGlobal("e", "energy")
    integrator.addComputePerDof("v", "a*v + (1-a)*f/m + b*sqrt(kT/m)*gaussian")
    integrator.addComputePerDof("x1", "x")
    integrator.addComputePerDof("x", "x+dt*v")
    integrator.addConstrainPositions()
    integrator.addComputePerDof("v", "(x-x1)/dt")
    return integrator


def time_steps(simulation, steps, repeats):
    samples = []
    for _ in range(repeats):
        simulation.context.getState(getEnergy=True)  # sync
        start = time.perf_counter()
        simulation.step(steps)
        simulation.context.getState(getEnergy=True)
        samples.append((time.perf_counter() - start) / steps * 1000.0)
    return statistics.median(samples), samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", default=None)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--only", default=None,
                        help="only time cases whose name contains this text")
    args = parser.parse_args()

    temperature = 300 * unit.kelvin
    forcefield, modeller = build_system(args.pdb)
    print("atoms: {}  platform: {}  steps x repeats: {} x {}".format(
        modeller.topology.getNumAtoms(), args.platform, args.steps,
        args.repeats))

    cases = [("reference Langevin (1 fused eval)", None)] + list(BOOSTS.items())
    if args.only:
        cases = [case for case in cases if args.only in case[0]]
    for name, create in cases:
        system = forcefield.createSystem(
            modeller.topology, nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds)
        if create is None:
            integrator = reference_integrator(temperature)
        else:
            integrator = create(system, temperature, **STAGE_ARGS)[2]
        simulation = app.Simulation(
            modeller.topology, system, integrator,
            openmm.Platform.getPlatformByName(args.platform))
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy(maxIterations=200)
        simulation.context.setVelocitiesToTemperature(temperature, 1)
        simulation.step(WARMUP_STEPS)
        median, samples = time_steps(simulation, args.steps, args.repeats)
        print("{:<36s} {:7.3f} ms/step  (runs: {})".format(
            name, median, ", ".join("{:.3f}".format(s) for s in samples)))
        del simulation


if __name__ == "__main__":
    main()
