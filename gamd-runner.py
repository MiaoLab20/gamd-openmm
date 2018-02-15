#!/usr/bin/python

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os
import traceback

from gamd import *


def createGamdSimulationFromAmberFiles(prmtopfile, inpcrdfile):
    prmtop = AmberPrmtopFile(prmtopfile)
    inpcrd = AmberInpcrdFile(inpcrdfile)
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = GamdIntegratorBoostTotalPotentialLowerBound(2.0 * femtoseconds, stage1end, stage2end, ntave,
                                                             6.0 * kilocalories_per_mole)
    simulation = Simulation(prmtop.topology, system, integrator)

    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    simulation.minimizeEnergy()

    return [simulation, integrator]

def createGamdSimulationFromPdbFile(pdbfile):
    pdb = PDBFile(pdbfile)
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                     nonbondedCutoff=1 * nanometer, constraints=HBonds)

    integrator = GamdIntegratorBoostTotalPotentialLowerBound(2.0 * femtoseconds, stage1end, stage2end, ntave,
                                                             6.0 * kilocalories_per_mole)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()

    return simulation



def createAmdSimulationFromAmberFiles(prmtopfile, inpcrdfile):
    prmtop = AmberPrmtopFile(prmtopfile)
    inpcrd = AmberInpcrdFile(inpcrdfile)
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = amd.AMDIntegrator(2.0*femtoseconds, 0, -1E99)
    simulation = Simulation(prmtop.topology, system, integrator)

    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    simulation.minimizeEnergy()

    return simulation


def createAmdSimulationFromPdbFile(pdbfile):
    pdb = PDBFile(pdbfile)
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                     nonbondedCutoff=1 * nanometer, constraints=HBonds)
    integrator = amd.AMDIntegrator(2.0 * femtoseconds, 0, -1E99)

    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()

    return simulation


def createLangevinSimulationFromAmberFiles(prmtopfile, inpcrdfile):
    prmtop = AmberPrmtopFile(prmtopfile)
    inpcrd = AmberInpcrdFile(inpcrdfile)
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    simulation = Simulation(prmtop.topology, system, integrator)

    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    simulation.minimizeEnergy()
    return simulation


ntave = 1000
stage1end = 2000
stage2end = 10000

os.makedirs("output/states/")
os.makedirs("output/checkpoints")

(simulation, integrator) = createGamdSimulationFromAmberFiles('dip.top', 'dip.crd')

simulation.saveState("output/states/validate.xml")
simulation.reporters.append(PDBReporter('output/output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, temperature=True,
                            potentialEnergy=True, totalEnergy=True, volume=True))

# Stage 1 - Conventional MD to gather statistics on Potential Values
simulation.step(stage1end)

# Debugging Version
# for step in range(stage1end):
#         simulation.step(1)
#         # step = count - 1
#         if step % ntave == (ntave - 1) or step % ntave == (ntave - 2) or step % ntave == 0 or step % ntave == 1:
#             integrator.printVariables()

# Stage 2
for step in range(stage1end + stage2end + 2):
    try:
        if step > 50:
            print("\nStep " + str(step) + " - prior to step:")
            integrator.printPositions()
        simulation.step(1)
        integrator.printVariables()
        simulation.saveState("output/states/" + str(step) + ".xml")
        simulation.saveCheckpoint("output/checkpoints/" + str(step) + ".bin")
    except Exception as err:
        print("\n\n-----------------------")
        print("Exception Information for")
        integrator.printVariables()
        print("\nException Message: " + str(err))
        traceback.print_exc(file=sys.stdout)
        integrator.printPositions()
        simulation.saveState("output/states/error.xml")
        simulation.saveCheckpoint("output/checkpoints/error.bin")
        break
