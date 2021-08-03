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
import datetime

from gamd.langevin.total_boost_integrators import LowerBoundIntegrator as TotalBoostLowerBoundIntegrator
from gamd.langevin.total_boost_integrators import UpperBoundIntegrator as TotalBoostUpperBoundIntegrator
from gamd.langevin.dihedral_boost_integrators import LowerBoundIntegrator as DihedralBoostLowerBoundIntegrator
from gamd.langevin.dihedral_boost_integrators import UpperBoundIntegrator as DihedralBoostUpperBoundIntegrator
from gamd.langevin.dual_boost_integrators import LowerBoundIntegrator as DualBoostLowerBoundIntegrator
from gamd.langevin.dual_boost_integrators import UpperBoundIntegrator as DualBoostUpperBoundIntegrator


def is_argument_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


def main():
    [boost_type, output_directory, device, platform] = handle_arguments()
    temperature = 298.15
    dt = 2.0 * femtoseconds
    ntcmdprep = 200000
    ntcmd = 1000000
    ntebprep = 200000
    nteb = 2000000
    nstlim = 18000000
    ntave = 25000
    number_of_steps_in_group = 100

    # This variable indicates the number of frames at the beginning of stage 5 (production) to ignore.
    #
    # starting_offset = 300
    starting_offset = 0

    start_date_time = datetime.datetime.now()
    print("Start Time: \t", start_date_time.strftime("%b-%d-%Y    %H:%M"))

    run_simulation(temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, boost_type, output_directory,
                   platform, device, number_of_steps_in_group, starting_offset)

    end_date_time = datetime.datetime.now()
    print("End Time: \t", end_date_time.strftime("%b-%d-%Y    %H:%M"))

    time_difference = end_date_time - start_date_time
    steps_per_second = nstlim / time_difference.seconds
    print("Execution rate for this run:  ", str(steps_per_second), " steps per second.")
    daily_execution_rate = steps_per_second * 3600 * 24 * dt
    print("Daily execution rate:         ", str(daily_execution_rate.value_in_unit(nanoseconds)), " ns per day.")
    production_starting_frame = ((ntcmd + nteb) / number_of_steps_in_group) + starting_offset
    run_post_simulation(temperature, output_directory, production_starting_frame)


#
#   NOTE:  Don't do this.  It moves the forces into separate groups, so that they don't get handled properly.
#
#    for i, force in enumerate(system.getForces()):
#        print(str(i) + "     " + force.__class__.__name__)
#        force.setForceGroup(i)
#        if force.__class__.__name__ == 'PeriodicTorsionForce':
#            group = i


def set_all_forces_to_group(system):
    group = 1
    for force in system.getForces():
        force.setForceGroup(group)
    return group


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


def create_lower_total_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    # The group is set, so that we can output the dihedral energy.  It doesn't impact calculations for total boost,
    # since we are utilizing the OpenMM provided variables with them not split out for total boost calculations.
    group = set_dihedral_group(system)
    return [group, TotalBoostLowerBoundIntegrator(dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                                       ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                                       ntave=ntave, temperature=temperature)]


def create_upper_total_boost_integrator(system, temperature, dt,ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    # The group is set, so that we can output the dihedral energy.  It doesn't impact calculations for total boost,
    # since we are utilizing the OpenMM provided variables with them not split out for total boost calculations.
    group = set_dihedral_group(system)
    return [group, TotalBoostUpperBoundIntegrator(dt=dt, ntcmdprep=ntcmdprep, ntcmd=ntcmd,
                                                                       ntebprep=ntebprep, nteb=nteb, nstlim=nstlim,
                                                                       ntave=ntave, sigma0=2.5 * kilocalories_per_mole, temperature=temperature)]


def create_lower_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave):
    group = set_all_forces_to_group(system)
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
    print("run-test.py boost-type [output-directory] [platform | device] [device]\n")
    print("\tboost-type:\t\tgamd-cmd-base|lower-total|upper-total|lower-dihedral|upper-dihedral|lower-dual|upper-dual\n")
    print("\toutput-directory:\tDirectory to output files. [default: output]\n")
    print("\tplatform:\t\tCUDA|OpenCL|CPU [defaults to OpenMM best guess for fastest]\n")
    print("\tdevice:\t\t\tUsed to specify the device index, when multiple GPUs exist on \n\t\t\t\tthe system. [defaults to CUDA, if platform not specified.]\n")


def handle_arguments():

    if len(sys.argv) == 2:
        output_directory = "output"
        boost_type = sys.argv[1]
        device = ""
        platform = ""
    elif len(sys.argv) == 3:
        boost_type = sys.argv[1]
        output_directory = sys.argv[2]
        device = ""
        platform = ""
    elif len(sys.argv) == 4:
        boost_type = sys.argv[1]
        output_directory = sys.argv[2]
        device = sys.argv[3]
        if is_argument_integer(device):
            platform = "CUDA"
        else:
            platform = device
            device = 0
    elif len(sys.argv) == 5:
        boost_type = sys.argv[1]
        output_directory = sys.argv[2]
        platform = sys.argv[3]
        device = sys.argv[4]
    else:
        usage()
        sys.exit(1)

    return [boost_type, output_directory, device, platform]


def run_simulation(unitless_temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, boost_type,
                   output_directory, platform_name, device, number_of_steps_in_group=100, reweighting_offset=0):
    coordinates_file = './data/md-4ns.rst7'
    prmtop_file = './data/dip.top'
    starttime = time.time()
    restarting = False
    restart_checkpoint_frequency = number_of_steps_in_group
    restart_checkpoint_filename = "gamd.backup"
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

    properties = {}
    if platform_name == "CUDA":
        platform = Platform.getPlatformByName(platform_name)
        properties['CudaPrecision'] = 'mixed'
        properties['DeviceIndex'] = device
        simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    elif platform_name == "OpenCL":
        platform = Platform.getPlatformByName(platform_name)
        properties['DeviceIndex'] = device
        simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    else:
        simulation = Simulation(prmtop.topology, system, integrator)
    
    """ # uncomment to view integrator computations
    for i in range(integrator.getNumComputations()):
        print(integrator.getComputationStep(i))
    
    print("exiting...")
    exit()
    """
    
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
        simulation.context.setVelocitiesToTemperature(unitless_temperature * kelvin)
        simulation.saveState(output_directory + "/states/initial-state.xml")
        simulation.reporters.append(DCDReporter(output_directory + '/output.dcd', number_of_steps_in_group))
        simulation.reporters.append(
            utils.ExpandedStateDataReporter(system, output_directory + '/state-data.log', number_of_steps_in_group,
                                            step=True,
                                            brokenOutForceEnergies=True, temperature=True, potentialEnergy=True,
                                            totalEnergy=True, volume=True))
        print("startup time: \t", time.time() - starttime)
        write_mode = "w"
        start_step = 1

    gamd_logger = GamdLogger(output_directory + "/gamd.log", write_mode, integrator, simulation)
    gamd_reweighting_logger = GamdLogger(output_directory + "/gamd-reweighting.log", write_mode, integrator, simulation)
    start_production_logging_step = ntcmd + nteb + (number_of_steps_in_group * reweighting_offset)
    start_production_logging_frame = int(start_production_logging_step / number_of_steps_in_group)

    with open(output_directory + "/" + "production-start-step.txt", "w") as prodstartstep_file:
        prodstartstep_file.write(str(start_production_logging_step))

    if not restarting:
        gamd_logger.write_header()
        gamd_reweighting_logger.write_header()

    print("Running: \t " + str(integrator.get_total_simulation_steps()) + " steps")
    for step in range(start_step, (integrator.get_total_simulation_steps() // number_of_steps_in_group) + 1):
        if step % restart_checkpoint_frequency // number_of_steps_in_group == 0:
            simulation.saveCheckpoint(restart_checkpoint_filename)

        # TEST
        #            if step == 250 and not restarting:
        #                print("sudden, unexpected interruption!")
        #                exit()
        # END TEST

        gamd_logger.mark_energies(group)
        if step >= start_production_logging_frame:
            gamd_reweighting_logger.mark_energies(group)

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
            if step >= start_production_logging_frame:
                gamd_reweighting_logger.write_to_gamd_log(step * number_of_steps_in_group)

            # print(integrator.get_current_state())
        except Exception as e:
            print("Failure on step " + str(step * number_of_steps_in_group))
            print(e)
            # print(integrator.get_current_state())
            gamd_logger.write_to_gamd_log(step)
            if step >= start_production_logging_frame:
                gamd_reweighting_logger.write_to_gamd_log(step * number_of_steps_in_group)

            sys.exit(2)

        # simulation.loadCheckpoint('/tmp/testcheckpoint')

        # debug_information = integrator.get_debugging_information()
        # getGlobalVariableNames(integrator)

        if step % integrator.ntave == 0:
            # if step % 1 == 0:

            simulation.saveState(output_directory + "/states/" + str(step * number_of_steps_in_group) + ".xml")
            simulation.saveCheckpoint(output_directory + "/checkpoints/" +
                                      str(step * number_of_steps_in_group) + ".bin")
            positions_filename = output_directory + '/positions/coordinates-' + \
                                 str(step * number_of_steps_in_group) + '.csv'
            integrator.create_positions_file(positions_filename)
            # pp = pprint.PrettyPrinter(indent=2)
            # pp.pprint(debug_information)


def create_graphics(execution_directory, command, temperature, output_filename):

    result = subprocess.run(["/bin/bash " + command + " " + " " + str(temperature)], capture_output=True, cwd=execution_directory, shell=True)

    with open(output_filename, "w") as output:
        output.write(result.stdout.decode('utf-8'))


def run_post_simulation(unitless_temperature, output_directory, starting_frame):
    with open(output_directory + "/"+ "temperature.dat", "w") as temperature_file:
        temperature_file.write(str(unitless_temperature))
    shutil.copy("tests/graphics/create-graphics.sh", output_directory + "/")
    write_out_cpptraj_command_files(output_directory, starting_frame)
    shutil.copytree("data", output_directory + "/data")
    create_graphics(output_directory + "/", "create-graphics.sh", str(unitless_temperature), 
                    output_directory + "/" + "graphics.log")


def write_out_cpptraj_command_files(output_directory, starting_frame):
    write_out_phi_cpptraj_command(output_directory, starting_frame)
    write_out_psi_cpptraj_command(output_directory, starting_frame)
    write_out_phi_psi_cpptraj_command(output_directory, starting_frame)


def write_out_phi_cpptraj_command(output_directory, starting_frame):
    with open(output_directory + "/" + "phi-dat-commands.cpptraj", "w") as dat_command_file:
        dat_command_file.write("trajin output.dcd " + str(int(starting_frame)) + "\n")
        dat_command_file.write("dihedral phi :1@C :2@N :2@CA :2@C out graphics/phi-cpptraj.dat" + "\n")
        dat_command_file.write("go" + "\n")


def write_out_psi_cpptraj_command(output_directory, starting_frame):
    with open(output_directory + "/" + "psi-dat-commands.cpptraj", "w") as dat_command_file:
        dat_command_file.write("trajin output.dcd " + str(int(starting_frame)) + "\n")
        dat_command_file.write("dihedral psi :2@N :2@CA :2@C :3@N out graphics/psi-cpptraj.dat" + "\n")
        dat_command_file.write("go" + "\n")


def write_out_phi_psi_cpptraj_command(output_directory, starting_frame):
    with open(output_directory + "/" + "phi-psi-commands.cpptraj", "w") as dat_command_file:
        dat_command_file.write("trajin output.dcd " + str(int(starting_frame)) + "\n")
        dat_command_file.write("dihedral phi :1@C :2@N :2@CA :2@C out graphics/phi-psi-cpptraj.dat" + "\n")
        dat_command_file.write("dihedral psi :2@N :2@CA :2@C :3@N out graphics/phi-psi-cpptraj.dat" + "\n")
        dat_command_file.write("go" + "\n")


if __name__ == "__main__":
    main()
