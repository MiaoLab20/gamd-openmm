from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from sys import exit
import os
import sys
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
    for i, force in enumerate(system.getForces()):
        force.setForceGroup(i)
        
    # group = 1
    # for force in system.getForces():
    #     print(force.__class__.__name__)
    #     force.setForceGroup(group)
    #     group += 1

#    integrator = LowerBoundIntegrator(2.0 * femtoseconds, 2000, 10000, 2000, 10000,
#                                                                     30000, 500)

    integrator = LowerBoundIntegrator()

    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    simulation.minimizeEnergy()

    simulation.saveState(output_directory + "/states/initial-state.xml")
    simulation.reporters.append(DCDReporter(output_directory + '/output.dcd', 1))
    simulation.reporters.append(utils.ExpandedStateDataReporter(system, output_directory + '/state-data.log', 1, step=True, brokenOutForceEnergies=True, temperature=True,
                                                  potentialEnergy=True, totalEnergy=True, volume=True))

    with open(output_directory + "/gamd.log", 'w') as gamdLog:
        gamdLog.write("total_nstep, Unboosted-Potential-Energy, Total-Force-Weight, Boost-Energy-Potential\n")
        
        for step in range(1, integrator.get_total_simulation_steps() + 1):
           try:
                # print(integrator.get_current_state())
               simulation.step(1)
               gamdLog.write(str(step) + "," + str(integrator.get_current_potential_energy()) + "," + str(integrator.get_force_scaling_factor()) + "," +
                             str(integrator.get_boost_potential()) + "\n")

                # print(integrator.get_current_state())
           except Exception as e:
               print("Failure on step " + str(step))
               print(e)
               print(integrator.get_current_state())
               gamdLog.write("Fail Step " + str(step) + "," + str(integrator.get_current_potential_energy()) + "," + str(integrator.get_force_scaling_factor()) + "," +
                              str(integrator.get_boost_potential()) + "\n")
               sys.exit(2)

           # debug_information = integrator.get_debugging_information()
           # getGlobalVariableNames(integrator)
        
           if step % integrator.ntave == 0:
            # if step % 1 == 0:

               simulation.saveState(output_directory + "/states/" + str(step) + ".xml")
               simulation.saveCheckpoint(output_directory + "/checkpoints/" + str(step) + ".bin")
               positions_filename = output_directory + '/positions/coordinates-' + str(step) + '.csv'
               integrator.create_positions_file(positions_filename)
               # pp = pprint.PrettyPrinter(indent=2)
               # pp.pprint(debug_information)

               
if __name__ == "__main__":
    main()
