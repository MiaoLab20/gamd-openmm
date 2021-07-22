#!/usr/bin/env python3

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
# from gamd.langevin.total_boost_integrators import UpperBoundIntegrator
# from gamd.langevin.dihedral_boost_integrators import import LowerBountIntegrator
# from gamd.langevin.dihedral_boost_integrators import import UpperBoundIntegrator
from gamd.stage_integrator import BoostType
from gamd.GamdLogger import GamdLogger
from gamd import utils as utils
import pprint
import shutil
import subprocess

from gamd.langevin.total_boost_integrators import LowerBoundIntegrator as TotalBoostLowerBoundIntegrator
from gamd.langevin.total_boost_integrators import UpperBoundIntegrator as TotalBoostUpperBoundIntegrator
from gamd.langevin.dihedral_boost_integrators import LowerBoundIntegrator as DihedralBoostLowerBoundIntegrator
from gamd.langevin.dihedral_boost_integrators import UpperBoundIntegrator as DihedralBoostUpperBoundIntegrator
from gamd.langevin.dual_boost_integrators import LowerBoundIntegrator as DualBoostLowerBoundIntegrator
from gamd.langevin.dual_boost_integrators import UpperBoundIntegrator as DualBoostUpperBoundIntegrator

def main():
    [boost_type, output_directory] = handle_arguments()
    temperature = 298.15
    dt = 2.0 * femtoseconds
    ntcmdprep = 200000
    ntcmd = 1000000
    ntebprep = 200000
    nteb = 1000000
    nstlim = 30000000
    ntave = 50000

    run_simulation(temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, boost_type, output_directory)
    run_post_simulation(temperature, output_directory)




#
#   NOTE:  Don't do this.  It moves the forces into separate groups, so that they don't get handled properly.
#
#    for i, force in enumerate(system.getForces()):
#        print(str(i) + "     " + force.__class__.__name__)
#        force.setForceGroup(i)
#        if force.__class__.__name__ == 'PeriodicTorsionForce':
#            group = i


def set_dihedral_group(system):
    group = 1
    for force in system.getForces():
        if force.__class__.__name__ == 'PeriodicTorsionForce':
            force.setForceGroup(group)
            break
    return group


def create_gamd_cmd_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep,
                                                                  nteb, nstlim, ntave):
    """
        This integrator is meant for use in generating a conventional MD baseline to compare against
        for the other integrators.

    :param system:
    :param temperature:
    :return:
    """
    group = set_dihedral_group(system)
    return [set_dihedral_group(system), DihedralBoostLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                                          ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                                          ntave=ntave, temperature=temperature,
                                                                          sigma0=0.0 * kilocalories_per_mole)]


def create_lower_total_boost_integrator(system, temperature, dt,ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    return [set_dihedral_group(system), TotalBoostLowerBoundIntegrator(dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                                       ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                                       ntave=ntave, temperature=temperature)]


def create_upper_total_boost_integrator(system, temperature, dt,ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    return [set_dihedral_group(system), TotalBoostUpperBoundIntegrator(dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                                       ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                                       ntave=ntave, temperature=temperature)]


def create_lower_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    group = set_dihedral_group(system)
    return [group, DihedralBoostLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                                     nteb=nteb, nstlim=nstlim, ntave=ntave, temperature=temperature)]


def create_upper_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    group = set_dihedral_group(system)
    return [group, DihedralBoostUpperBoundIntegrator(group,  dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                                     nteb=nteb, nstlim=nstlim, ntave=ntave, temperature=temperature)]
    
def create_lower_dual_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    group = set_dihedral_group(system)
    return [group, DualBoostLowerBoundIntegrator(group, dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                                     nteb=nteb, nstlim=nstlim, ntave=ntave, temperature=temperature)]


def create_upper_dual_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    group = set_dihedral_group(system)
    return [group, DualBoostUpperBoundIntegrator(group,  dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd, ntebprep=ntebprep,
                                                     nteb=nteb, nstlim=nstlim, ntave=ntave, temperature=temperature)]


def create_output_directories(directories):
    for dir in directories:
        os.makedirs(dir, 0o755)


def usage():
    print("run-test.py boost-type [output-directory]\n")
    print("\tboost-type:\t\tgamd-cmd-base|lower-total|upper-total|lower-dihedral|upper-dihedral|lower-dual|upper-dual\n")
    print("\toutput-directory:\tDirectory to output files. [default: output]\n")


def handle_arguments():

    if len(sys.argv) == 2:
        output_directory = "output"
        boost_type = sys.argv[1]
    elif len(sys.argv) == 3:
        boost_type = sys.argv[1]
        output_directory = sys.argv[2]
    else:
        usage()
        sys.exit(1)

    return [boost_type, output_directory]


def run_simulation(unitless_temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, boost_type, output_directory):
    coordinates_file = './data/md-4ns.rst7'
    prmtop_file = './data/dip.top'
    starttime = time.time()
    restarting = False
    restart_checkpoint_frequency = 1000
    restart_checkpoint_filename = "gamd.backup"
    number_of_steps_in_group = 1000
    temperature = unitless_temperature * kelvin
    prmtop = AmberPrmtopFile(prmtop_file)
    inpcrd = AmberInpcrdFile(coordinates_file)
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.8 * nanometer, constraints = HBonds)
    # dihedral_boost = True
    # (simulation, integrator) = createGamdSimulationFromAmberFiles(prmtop_file, coordinates_file, dihedral_boost=dihedral_boost)


    if boost_type == "gamd-cmd-base":
        [group, integrator] = create_gamd_cmd_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb,
                                                         nstlim, ntave)
    elif boost_type == "lower-total":
        [group, integrator] = create_lower_total_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep,
                                                                  nteb, nstlim, ntave)
    elif boost_type == "upper-total":
        [group, integrator] = create_upper_total_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep,
                                                                  nteb, nstlim, ntave)
    elif boost_type == "lower-dihedral":
        [group, integrator] = create_lower_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep,
                                                                  nteb, nstlim, ntave)
    elif boost_type == "upper-dihedral":
        [group, integrator] = create_upper_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep,
                                                                  nteb, nstlim, ntave)
    elif boost_type == "lower-dual":
        [group, integrator] = create_lower_dual_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep,
                                                                  nteb, nstlim, ntave)
    elif boost_type == "upper-dual":
        [group, integrator] = create_upper_dual_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep,
                                                                  nteb, nstlim, ntave)
        
    else:
        usage()
        sys.exit(1)

    if not restarting:
        create_output_directories([output_directory, output_directory + "/states/", output_directory + "/positions/",
                                   output_directory + "/pdb/", output_directory + "/checkpoints"])


    simulation = Simulation(prmtop.topology, system, integrator)
    if restarting:
        simulation.loadCheckpoint(restart_checkpoint_filename)
        write_mode = "a"
        start_step = int(integrator.getGlobalVariableByName("stepCount") // number_of_steps_in_group)
        print("restarting from saved checkpoint:", restart_checkpoint_filename,
              "at step:", start_step)
    else:
        simulation.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        simulation.minimizeEnergy()
        #
        # Should we be setting this to the temperature?
        #
        simulation.context.setVelocitiesToTemperature(298.15 * kelvin)
        simulation.saveState(output_directory + "/states/initial-state.xml")
        simulation.reporters.append(DCDReporter(output_directory + '/output.dcd', number_of_steps_in_group))
        simulation.reporters.append(
            utils.ExpandedStateDataReporter(system, output_directory + '/state-data.log', number_of_steps_in_group,
                                            step=True,
                                            brokenOutForceEnergies=True, temperature=True, potentialEnergy=True,
                                            totalEnergy=True, volume=True))
        print("startup time:", time.time() - starttime)
        write_mode = "w"
        start_step = 1

    gamd_logger = GamdLogger(output_directory + "/gamd.log", write_mode, integrator, simulation)

    if not restarting:
        gamd_logger.write_header()
    print("Running " + str(integrator.get_total_simulation_steps()) + " steps")
    for step in range(start_step, (integrator.get_total_simulation_steps() // number_of_steps_in_group) + 1):
        if step % restart_checkpoint_frequency // number_of_steps_in_group == 0:
            simulation.saveCheckpoint(restart_checkpoint_filename)

        # TEST
        #            if step == 250 and not restarting:
        #                print("sudden, unexpected interruption!")
        #                exit()
        # END TEST

        gamd_logger.mark_energies(group)
        try:
            # print(integrator.get_current_state())

            #
            #  NOTE:  We need to save off the starting total and dihedral potential energies, since we
            #         end up modifying them by updating the state of the particles.  This allows us to write
            #         out the values as they were at the beginning of the step for what all of the calculations
            #         for boosting were based on.
            #

            simulation.step(number_of_steps_in_group)
            gamd_logger.write_to_gamd_log(step * number_of_steps_in_group)

            # print(integrator.get_current_state())
        except Exception as e:
            print("Failure on step " + str(step * number_of_steps_in_group))
            print(e)
            # print(integrator.get_current_state())
            gamd_logger.write_to_gamd_log(step)
            sys.exit(2)

        # simulation.loadCheckpoint('/tmp/testcheckpoint')

        # debug_information = integrator.get_debugging_information()
        # getGlobalVariableNames(integrator)

        if step % integrator.ntave == 0:
            # if step % 1 == 0:

            simulation.saveState(output_directory + "/states/" + str(step * number_of_steps_in_group) + ".xml")
            simulation.saveCheckpoint(output_directory + "/checkpoints/" + \
                                      str(step * number_of_steps_in_group) + ".bin")
            positions_filename = output_directory + '/positions/coordinates-' + \
                                 str(step * number_of_steps_in_group) + '.csv'
            integrator.create_positions_file(positions_filename)
            # pp = pprint.PrettyPrinter(indent=2)
            # pp.pprint(debug_information)


def create_graphics(execution_directory, command, temperature, output_filename):
    result = subprocess.run([command, temperature], cwd=execution_directory, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    with open(output_filename, "w") as output:
        output.write(result.stdout.decode('utf-8'))


def run_post_simulation(unitless_temperature, output_directory):
    shutil.copy("tests/graphics/create-graphics.sh", output_directory + "/")
    shutil.copy("tests/graphics/phi-dat-commands.cpptraj", output_directory + "/")
    shutil.copy("tests/graphics/psi-dat-commands.cpptraj", output_directory + "/")
    shutil.copy("tests/graphics/phi-psi-commands.cpptraj", output_directory + "/")
    create_graphics(output_directory, "create-graphics.sh", str(unitless_temperature), "graphics.log")


if __name__ == "__main__":
    main()
