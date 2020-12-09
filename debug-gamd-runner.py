from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from sys import exit
import os
import sys
import time
import traceback
# from gamd.langevin.total_boost_integrators import LowerBoundIntegrator
from gamd.langevin.total_boost_integrators import UpperBoundIntegrator
# from gamd.langevin.dihedral_boost_integrators import DihedralLowerBoundIntegrator
# from gamd.langevin.dihedral_boost_integrators import UpperBoundIntegrator
from gamd.stage_integrator import BoostType
from gamd import utils as utils
import pprint


def create_output_directories(directories):
    for dir in directories:
        os.makedirs(dir, 0o755)


def get_global_variable_names(integrator):
    for index in range(0, integrator.getNumGlobalVariables()):
        print(integrator.getGlobalVariableName(index))


def print_global_variables(integrator):
    for index in range(0, integrator.getNumGlobalVariables()):
        name = integrator.getGlobalVariableName(index)
        value = integrator.getGlobalVariableByName(name)
        print(name + ":  " + str(value))


def print_global_variables_headers(integrator, output):
    number_of_globals = integrator.getNumGlobalVariables()
    for index in range(0, number_of_globals):
        name = integrator.getGlobalVariableName(index)
        output.write(str(name))
        if index < number_of_globals:
            output.write(",")
    output.write("\n")


def print_global_variables_values(integrator, output):
    number_of_globals = integrator.getNumGlobalVariables()
    for index in range(0, number_of_globals):
        name = integrator.getGlobalVariableName(index)
        value = integrator.getGlobalVariableByName(name)
        output.write(str(value))
        if index < number_of_globals:
            output.write(",")
    output.write("\n")


def write_to_gamd_log(gamd_log, step, group, simulation, integrator):
    state = simulation.context.getState(getEnergy=True)
    dihedral_state = simulation.context.getState(getEnergy=True, groups={group})
    force_scaling_factors = integrator.get_force_scaling_factors()
    boost_potentials = integrator.get_boost_potentials()

    total_potential_energy = str(state.getPotentialEnergy() / (kilojoules_per_mole * 4.184))
    dihedral_energy = str(dihedral_state.getPotentialEnergy() / (kilojoules_per_mole * 4.184))
    total_force_scaling_factor = str(force_scaling_factors[BoostType.TOTAL.value + "ForceScalingFactor"])
    dihedral_force_scaling_factor = str(force_scaling_factors[BoostType.DIHEDRAL.value + "ForceScalingFactor"])
    total_boost_potential = str(boost_potentials[BoostType.TOTAL.value + "BoostPotential"] / 4.184)
    dihedral_boost_potential = str(boost_potentials[BoostType.DIHEDRAL.value + "BoostPotential"] / 4.184)

    gamd_log.write("\t" + str(1) + "\t" + str(step * 1) + "\t" +
                   total_potential_energy + "\t" +
                   dihedral_energy + "\t" +
                   total_force_scaling_factor + "\t" +
                   dihedral_force_scaling_factor + "\t" +
                   total_boost_potential + "\t" +
                   dihedral_boost_potential + "\n")


def main():
    starttime = time.time()
    if len(sys.argv) == 1:
        output_directory = "output"
    else:
        output_directory = sys.argv[1]

    coordinates_file = './data/md-4ns.rst7'
    prmtop_file = './data/dip.top'
    
    restarting = False
    restart_checkpoint_frequency = 1000
    restart_checkpoint_filename = "gamd.backup"
    
    if not restarting:
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

    #
    #   NOTE:  Don't do this.  It moves the forces into seperate groups, so that they don't get handled properly.
    #
    #    for i, force in enumerate(system.getForces()):
    #        print(str(i) + "     " + force.__class__.__name__)
    #        force.setForceGroup(i)
    #        if force.__class__.__name__ == 'PeriodicTorsionForce':
    #            group = i

    # Total Boost
    # integrator = LowerBoundIntegrator(dt=2.0 * femtoseconds, ntcmdprep=2000, ntcmd=10000, ntebprep=2000,
    #                                  nteb=10000, nstlim=30000, ntave=500)
    integrator = UpperBoundIntegrator(dt=2.0 * femtoseconds, ntcmdprep=2000, ntcmd=10000, ntebprep=2000,
                                      nteb=10000, nstlim=30000, ntave=500)

    # Dihedral Boost
    #integrator = DihedralLowerBoundIntegrator(group, dt=2.0 * femtoseconds, ntcmdprep=2000, ntcmd=10000, ntebprep=2000,
    #                                  nteb=10000, nstlim=30000, ntave=500)
    # integrator = UpperBoundIntegrator(group, dt=2.0 * femtoseconds, ntcmdprep=2000, ntcmd=10000, ntebprep=2000,
    #                                   nteb=10000, nstlim=30000, ntave=500)


    simulation = Simulation(prmtop.topology, system, integrator)
    if restarting:
        simulation.loadCheckpoint(restart_checkpoint_filename)
        write_mode = "a"
        start_step = int(integrator.getGlobalVariableByName("stepCount") // 100)
        print("restarting from saved checkpoint:", restart_checkpoint_filename,
              "at step:", start_step)
    else:
        simulation.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(298.15*kelvin)
        simulation.saveState(output_directory + "/states/initial-state.xml")
        simulation.reporters.append(DCDReporter(output_directory + '/output.dcd', 100))
        simulation.reporters.append(utils.ExpandedStateDataReporter(system, output_directory + '/state-data.log', 1, step=True, brokenOutForceEnergies=True, temperature=True,
                                                      potentialEnergy=True, totalEnergy=True, volume=True))
        print("startup time:", time.time() - starttime)
        write_mode = "w"
        start_step = 1
    with open(output_directory + "/gamd.log", write_mode) as gamdLog:
        with open(output_directory + "/debug.log", write_mode) as debugLog:
            print_global_variables_headers(integrator, debugLog)

            if not restarting:
                gamdLog.write("# Gaussian accelerated Molecular Dynamics log file\n")
                gamdLog.write("# All energy terms are stored in unit of kcal/mol\n")
                gamdLog.write("# ntwx,total_nstep,Unboosted-Potential-Energy,Unboosted-Dihedral-Energy,Total-Force-Weight,Dihedral-Force-Weight,Boost-Energy-Potential,Boost-Energy-Dihedral\n")

            for step in range(start_step, (integrator.get_total_simulation_steps() // 1) + 1):
                if step % restart_checkpoint_frequency // 100 == 0:
                    simulation.saveCheckpoint(restart_checkpoint_filename)

                # TEST
    #            if step == 250 and not restarting:
    #                print("sudden, unexpected interruption!")
    #                exit()
                # END TEST

                try:
                    # print(integrator.get_current_state())
                    simulation.step(1)
                    print_global_variables_values(integrator, debugLog)
                    write_to_gamd_log(gamdLog, step, group, simulation, integrator)

                    # print(integrator.get_current_state())
                except Exception as e:
                    print("Failure on step " + str(step*1))
                    print(e)
                    print("BLOWING UP!!!")
                    print_global_variables(integrator)
                    print_global_variables_values(integrator, debugLog)
                    # print(integrator.get_current_state())
                    write_to_gamd_log(gamdLog, step, group, simulation, integrator)
                    sys.exit(2)


                #simulation.loadCheckpoint('/tmp/testcheckpoint')

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
