from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from sys import exit
import os
import sys
import time
import traceback
from gamd.langevin.total_boost_integrators import LowerBoundIntegrator
from gamd import utils as utils
import pprint


def create_output_directories(directories):
    for dir in directories:
        os.makedirs(dir, 0o755)


def getGlobalVariableNames(integrator):
    for index in range(0, integrator.getNumGlobalVariables()):
        print(integrator.getGlobalVariableName(index))


def main():
    starttime = time.time()
    if len(sys.argv) == 1:
        output_directory = "output"
    else:
        output_directory = sys.argv[1]

    coordinates_file = './data/md-4ns.rst7'
    prmtop_file = './data/dip.top'

    create_output_directories([output_directory, output_directory + "/states/", output_directory + "/positions/",
                               output_directory + "/pdb/", output_directory + "/checkpoints"])

    # dihedral_boost = True

    # (simulation, integrator) = createGamdSimulationFromAmberFiles(prmtop_file, coordinates_file, dihedral_boost=dihedral_boost)

    prmtop = AmberPrmtopFile(prmtop_file)
    inpcrd = AmberInpcrdFile(coordinates_file)
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=HBonds)
    # for i, force in enumerate(system.getForces()):
    #    force.setForceGroup(i)
        
    group = 1
    for force in system.getForces():
    #     print(force.__class__.__name__)
        if force.__class__.__name__ == 'PeriodicTorsionForce':
            force.setForceGroup(group)
            break
    #     group += 1

#    integrator = LowerBoundIntegrator(2.0 * femtoseconds, 2000, 10000, 2000, 10000,
#                                                                     30000, 500)

    integrator = LowerBoundIntegrator()
    #integrator = UpperBoundIntegrator()

    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(298.15*kelvin)
    simulation.saveState(output_directory + "/states/initial-state.xml")
    simulation.reporters.append(DCDReporter(output_directory + '/output.dcd', 100))
    simulation.reporters.append(utils.ExpandedStateDataReporter(system, output_directory + '/state-data.log', 100, step=True, brokenOutForceEnergies=True, temperature=True,
                                                  potentialEnergy=True, totalEnergy=True, volume=True))
    print("startup time:", time.time() - starttime)
    with open(output_directory + "/gamd.log", 'w') as gamdLog:
        gamdLog.write("# Gaussian accelerated Molecular Dynamics log file\n")
        gamdLog.write("# All energy terms are stored in unit of kcal/mol\n")
        gamdLog.write("# ntwx,total_nstep,Unboosted-Potential-Energy,Unboosted-Dihedral-Energy,Total-Force-Weight,Dihedral-Force-Weight,Boost-Energy-Potential,Boost-Energy-Dihedral\n")

        for step in range(1, (integrator.get_total_simulation_steps() // 100) + 1):
           try:
               # print(integrator.get_current_state())
               simulation.step(100)
               state = simulation.context.getState(getEnergy=True, groups={group})
               gamdLog.write("\t" + str(100) + "\t" + str(step*100) + "\t" + str(integrator.get_current_potential_energy()/4.184) + "\t" + str(state.getPotentialEnergy()/(kilojoules_per_mole*4.184)) + "\t" +  str(integrator.get_total_force_scaling_factor()) + "\t" + str(integrator.get_dihedral_force_scaling_factor()) + "\t" +
                             str(integrator.get_boost_potential()/4.184) + "\t" +  str(integrator.get_dihedral_boost()/4.184) + "\n")

                # print(integrator.get_current_state())
           except Exception as e:
               print("Failure on step " + str(step*100))
               print(e)
               print(integrator.get_current_state())
               state = simulation.context.getState(getEnergy=True, groups={group})
               gamdLog.write("Fail Step: " + str(step*100) + "\t" + str(integrator.get_current_potential_energy()/4.184) + "\t" + str(state.getPotentialEnergy()/(4.184*kilojoules_per_mole)) + "\t" + str(integrator.get_total_force_scaling_factor()) + "\t" + str(integrator.get_dihedral_force_scaling_factor()) + "\t" +
                              str(integrator.get_boost_potential()/4.184) + "\t" + str(integrator.get_dihedral_boost()/4.184) + "\n")
               sys.exit(2)

           # debug_information = integrator.get_debugging_information()
           # getGlobalVariableNames(integrator)
        
           if step % integrator.ntave == 0:
            # if step % 1 == 0:

               simulation.saveState(output_directory + "/states/" + str(step*100) + ".xml")
               simulation.saveCheckpoint(output_directory + "/checkpoints/" + str(step*100) + ".bin")
               positions_filename = output_directory + '/positions/coordinates-' + str(step*100) + '.csv'
               integrator.create_positions_file(positions_filename)
               # pp = pprint.PrettyPrinter(indent=2)
               # pp.pprint(debug_information)

               
if __name__ == "__main__":
    main()
