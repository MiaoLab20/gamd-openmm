#!/usr/bin/python3

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from sys import exit
import os
import traceback

from gamd import *


def createPositionsFile(integrator, thefile):
    positions = integrator.get_coordinates()
    f = open(thefile, 'w')
    f.write("particle, x, y, z\n")
    for i in range(len(positions)):
       f.write(str(i) + ", " + str(positions[i][0]) + ", " + str(positions[i][1])+ ", " + str(positions[i][2]) +"\n" )
    f.close()

def createGamdLog(gamdLog, filename):
    with  open(filename, 'w') as f:
        keys = gamdLog[0].keys()
        for header in keys[:-1]:
            f.write(header + ", ")
        f.write(keys[-1] + "\n")
        for entry in gamdLog:
            for header in keys[:-1]:
                f.write(str(entry[header]) + ", ")
            f.write(str(entry[keys[-1]]) + "\n")


def createGamdSimulationFromAmberFiles(prmtopfile, inpcrdfile, lowerBound=True, dihedral_boost=False):
    prmtop = AmberPrmtopFile(prmtopfile)
    inpcrd = AmberInpcrdFile(inpcrdfile)
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    if dihedral_boost:
        # Move torsion terms to force group 1
        group = 1
        for force in system.getForces():
            if force.__class__.__name__ != 'PeriodicTorsionForce':
                force.setForceGroup(group)
                break
        
        if lowerBound:
            integrator = GamdGroupBoostPotentialIntegratorLowerBound(2.0 * femtoseconds, group, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole)
        else:
            integrator = GamdGroupBoostPotentialIntegratorUpperBound(2.0 * femtoseconds, group, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole)
    else:
        if lowerBound:
            integrator = GamdTotalBoostPotentialIntegratorLowerBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole)
        else:
            integrator = GamdTotalBoostPotentialIntegratorUpperBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole)

    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    simulation.minimizeEnergy()

    return [simulation, integrator]


def createGamdSimulationFromPdbFile(pdbfile, prmtopfile, lowerBound=True, dihedral_boost=False):
    pdb = PDBFile(pdbfile)
    prmtop = AmberPrmtopFile(prmtopfile)
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1 * nanometer, constraints=HBonds)
    if dihedral_boost:
        group = 1
        for force in system.getForces():
            if force.__class__.__name__ != 'PeriodicTorsionForce':
                force.setForceGroup(group)
                break
        
        if lowerBound:
            integrator = GamdGroupBoostPotentialIntegratorLowerBound(2.0 * femtoseconds, group, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole)
        else:
            integrator = GamdGroupBoostPotentialIntegratorUpperBound(2.0 * femtoseconds, group, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole)
    else:
        if lowerBound:
            integrator = GamdTotalBoostPotentialIntegratorLowerBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole)
        else:
            integrator = GamdTotalBoostPotentialIntegratorUpperBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole)
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

    simulation.minimizeEnergy()

    return [simulation, integrator]


def create_output_directories(directories):
    for dir in directories:
        os.makedirs(dir, 0o755)


def upperBoundFunctionWrapper(a_function):
    def return_function(prmtopfile, coordinatesfile):
        a_function(prmtopfile, coordinatesfile, False)
    return return_function
  
def lowerBoundGroupFunctionWrapper(a_function):
    def return_function(prmtopfile, coordinatesfile):
        a_function(prmtopfile, coordinatesfile, False, True)
    return return_function
  
def upperBoundGroupFunctionWrapper(a_function):
    def return_function(prmtopfile, coordinatesfile):
        a_function(prmtopfile, coordinatesfile, False, True)
    return return_function


def main():

    output_directory = "output"
    coordinates_file = './data/md-4ns.rst7'
    prmtop_file =  './data/dip.top'

    total_number_of_steps = number_of_steps_in_stage_1 + number_of_steps_in_stage_2 + number_of_steps_in_stage_3

    function_dictionary = {'total-lb': {'amber' : createGamdSimulationFromAmberFiles,
                                        'pdb': createGamdSimulationFromPdbFile},
                            'total-ub': {'amber' : upperBoundFunctionWrapper(createGamdSimulationFromAmberFiles),
                                         'pdb': upperBoundFunctionWrapper(createGamdSimulationFromPdbFile) },
                           'group-lb': {'amber' : lowerBoundGroupFunctionWrapper(createGamdSimulationFromAmberFiles),
                                        'pdb': lowerBoundGroupFunctionWrapper(createGamdSimulationFromPdbFile)},
                            'group-ub': {'amber' : upperBoundGroupFunctionWrapper(createGamdSimulationFromAmberFiles),
                                         'pdb': upperBoundGroupFunctionWrapper(createGamdSimulationFromPdbFile) }}


    create_output_directories([output_directory, output_directory + "/states/", output_directory + "/positions/",
                               output_directory + "/pdb/", output_directory + "/checkpoints"])

    dihedral_boost = True
    (simulation, integrator) = createGamdSimulationFromAmberFiles(prmtop_file, coordinates_file, dihedral_boost=dihedral_boost)


    simulation.saveState(output_directory + "/states/initial-state.xml")
    simulation.reporters.append(PDBReporter(output_directory + '/output.pdb', 10000))
    simulation.reporters.append(StateDataReporter(stdout, 10000, step=True, temperature=True,
                                                  potentialEnergy=True, totalEnergy=True, volume=True))
    gamdLog = []
    for step in range(total_number_of_steps):
        simulation.step(1)
        if step % ntave == 0:
            simulation.saveState(output_directory + "/states/" + str(step) + ".xml")
            simulation.saveCheckpoint(output_directory + "/checkpoints/" + str(step) + ".bin")
            gamdLog.append({'total_nstep': step,
                            'Unboosted-Potential-Energy': integrator.get_current_potential_energy(),
                            'Total-Force-Weight': integrator.get_total_force_scaling_factor(),
                            'Boost-Energy-Potential': integrator.get_boost_potential() })
            createPositionsFile(integrator, output_directory + '/positions/coordinates-' + str(step) + '.csv')

    createGamdLog(gamdLog, output_directory + "/gamd.log")

ntave = 1000
number_of_steps_in_stage_1 = 10000
number_of_steps_in_stage_2 = 100000
number_of_steps_in_stage_3 = 100000

if __name__ == "__main__":
    main()
