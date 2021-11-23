#!/usr/bin/env python3

import datetime
import shutil
import subprocess
import time

import parmed
from simtk.openmm import *
from simtk.openmm.app import *

from gamd import utils as utils
from gamd.DebugLogger import DebugLogger
from gamd.GamdLogger import GamdLogger
from gamd.integrator_factory import *


class RunningRates:

    def __init__(self, number_of_simulation_steps: int, save_rate: int, reporting_rate: int,
                 debugging_enabled: bool = False, debugging_step_function=None) -> None:
        """Constructor

            The save_rate determines that rate at which to write out
            DCD/PDB file, checkpoint, coordinates, and write to the GaMD log files.
            The value of save_rate should evenly divide the number of simulation steps.

            The reporting_rate determines the rate at which to write out to
            the state-data.log and the debugging files, if enabled.  The value of the
            reporting rate should evenly divide the number of simulation steps.

            The save_rate or reporting_rate should be a multiple of the other
            value, if debugging is enabled.

            The save_rate, reporting_rate, and debugging flag determine the
            batch_run_rate, which is just the number of steps to have
            OpenMM execute at one time in the simulation.

        Parameters

        :param number_of_simulation_steps: (int) Same as nstlim, the total number of steps in the simulation.

        :param save_rate: (int) the number of steps before performing a save.

        :param reporting_rate: (int) the number of steps before printing to reports.

        :param debugging_enabled: (boolean) indicates whether debugging is enabled.  Defaults: False

        """
        if ((number_of_simulation_steps % save_rate) != 0 and
                (number_of_simulation_steps % reporting_rate) != 0):
            raise ValueError("RunningRates:  The save_rate and reporting_rate should evenly "
                             "divide into the number of steps in the simulation.")

        if debugging_enabled and not (((save_rate % reporting_rate) == 0) or
                                      ((reporting_rate % save_rate) == 0)):
            raise ValueError("RunningRates:  When debugging is enabled, save_rate/reporting_rate "
                             "must be a multiple of the other value.")

        self.debugging_step_function = debugging_step_function
        if callable(debugging_step_function):
            self.custom_debugging_function = True
        else:
            self.custom_debugging_function = False

        self.save_rate = save_rate
        self.reporting_rate = reporting_rate
        self.number_of_simulation_steps = number_of_simulation_steps
        self.debugging_enabled = debugging_enabled

        if debugging_enabled:
            if save_rate <= reporting_rate:
                self.batch_run_rate = save_rate
            else:
                self.batch_run_rate = reporting_rate
        else:
            self.batch_run_rate = self.save_rate

    def get_save_rate(self):
        return self.save_rate

    def get_reporting_rate(self):
        return self.reporting_rate

    def is_save_step(self, step):
        return (step % self.save_rate) == 0

    def is_reporting_step(self, step):
        return (step % self.reporting_rate) == 0

    def get_batch_run_rate(self):
        return self.batch_run_rate

    def is_debugging_step(self, step):
        if self.custom_debugging_function:
            result = self.debugging_step_function(step)
        else:
            result = (step % self.reporting_rate) == 0
        return result

    def get_batch_run_range(self):
        return range(1, self.number_of_simulation_steps // self.get_batch_run_rate() + 1)

    def get_restart_step(self, integrator):
        start_step = int(integrator.getGlobalVariableByName("stepCount") // self.get_batch_run_rate())
        return start_step

    def get_restart_batch_run_range(self, integrator):
        return range(self.get_restart_step(integrator),
                     self.number_of_simulation_steps // self.get_batch_run_rate() + 1)

    def get_step_from_frame(self, frame):
        return frame * self.get_batch_run_rate()


def is_argument_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


def main():
    [boost_type, output_directory, device, platform, debug, quick] = handle_arguments()
    temperature = 300
    dt = 2.0 * femtoseconds

    if quick:
        ntcmdprep = 2000
        ntcmd = 10000
        ntebprep = 2000
        nteb = 20000
        nstlim = 60000
        ntave = 250
        frame_size = 50
    elif boost_type == "gamd-cmd-base":
        ntcmdprep = 200000
        ntcmd = 1000000
        ntebprep = 200000
        nteb = 2000000
        nstlim = 500000000
        ntave = 25000
        frame_size = 100
    else:
        ntcmdprep = 200000
        ntcmd = 1000000
        ntebprep = 200000
        nteb = 2000000
        nstlim = 18000000
        ntave = 25000
        frame_size = 10

    if debug:
        running_rates = RunningRates(nstlim, frame_size, 1, True)
    else:
        running_rates = RunningRates(nstlim, frame_size, frame_size, False)


    # This variable indicates the number of frames at the beginning of stage 5 (production) to ignore.
    #
    # starting_offset = 300
    starting_offset = 0

    start_date_time = datetime.datetime.now()
    print("Start Time: \t", start_date_time.strftime("%b-%d-%Y    %H:%M"))

    start_date_time = run_simulation(temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, boost_type,
                                     output_directory, platform, device, running_rates, starting_offset,
                                     debug)


    output_starting_parameters(output_directory, temperature, dt, ntcmdprep, ntcmd, ntebprep,
                               nteb, nstlim, ntave, frame_size)


    print("Start Time: \t", start_date_time.strftime("%b-%d-%Y    %H:%M"))
    end_date_time = datetime.datetime.now()
    print("End Time: \t", end_date_time.strftime("%b-%d-%Y    %H:%M"))

    time_difference = end_date_time - start_date_time
    steps_per_second = nstlim / time_difference.seconds
    print("Execution rate for this run:  ", str(steps_per_second), " steps per second.")
    daily_execution_rate = steps_per_second * 3600 * 24 * dt
    print("Daily execution rate:         ", str(daily_execution_rate.value_in_unit(nanoseconds)), " ns per day.")
    production_starting_frame = ((ntcmd + nteb) / frame_size) + starting_offset
    run_post_simulation(temperature, output_directory, production_starting_frame)

#
#   NOTE:  Don't do this.  It moves the forces into separate groups, so that they don't get handled properly.
#
#    for i, force in enumerate(system.getForces()):
#        print(str(i) + "     " + force.__class__.__name__)
#        force.setForceGroup(i)
#        if force.__class__.__name__ == 'PeriodicTorsionForce':
#            group = i


def output_starting_parameters(output_directory, temperature, dt, ntcmdprep, ntcmd,
                               ntebprep, nteb, nstlim, ntave, frame_size):
    filename = os.path.join(output_directory, "starting-parameters.txt")
    with open(filename, "w") as output:
        output.write("temperature={0}\n".format(temperature))
        output.write("dt={0}\n".format(dt))
        output.write("ntcmdprep={0}\n".format(ntcmdprep))
        output.write("ntcmd={0}\n".format(ntcmd))
        output.write("ntebprep={0}\n".format(ntebprep))
        output.write("nteb={0}\n".format(nteb))
        output.write("nstlim={0}\n".format(nstlim))
        output.write("ntave={0}\n".format(ntave))
        output.write("frame_size={0}\n".format(frame_size))




def print_integration_algorithms(filename):
    #pdb_filename = './data/rna_ds_separated25A_wet.pdb'
    pdb_filename = './data/rna_ds_separated25A_dry.pdb'
    pdb = PDBFile(pdb_filename)
    pdb_parmed_for_box = parmed.load_file(pdb_filename)
    # Use for explicit solvent
    #forcefield = ForceField('amoeba2018.xml')
    # Use for implicit solvent
    forcefield = ForceField('amoeba2018.xml', 'amoeba2018_gk.xml') # for implicit solvent
    # Use for explicit solvent
    #system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=HBonds)
    # Use for implicit solvent
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    temperature = unitless_temperature
    dt = 2.0 * femtoseconds
    ntcmdprep = 200000
    ntcmd = 1000000
    ntebprep = 200000
    nteb = 2000000
    nstlim = 18000000
    ntave = 25000

    [group, group, integrator] = create_lower_total_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                                     ntebprep,
                                                                     nteb, nstlim, ntave)
    DebugLogger.write_integration_algorithm_to_file(filename, integrator)

    [group, group, integrator] = create_upper_total_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                                     ntebprep,
                                                                     nteb, nstlim, ntave)
    DebugLogger.write_integration_algorithm_to_file(filename, integrator)

    [group, group, integrator] = create_lower_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                                        ntebprep,
                                                                        nteb, nstlim, ntave)
    DebugLogger.write_integration_algorithm_to_file(filename, integrator)

    [group, group, integrator] = create_upper_dihedral_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                                        ntebprep,
                                                                        nteb, nstlim, ntave)
    DebugLogger.write_integration_algorithm_to_file(filename, integrator)

    [group, group, integrator] = create_lower_dual_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep,
                                                                    nteb, nstlim, ntave)
    DebugLogger.write_integration_algorithm_to_file(filename, integrator)

    [group, group, integrator] = create_upper_dual_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd, ntebprep,
                                                                    nteb, nstlim, ntave)
    DebugLogger.write_integration_algorithm_to_file(filename, integrator)

    [group, group, integrator] = create_lower_non_bonded_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                                          ntebprep,
                                                                          nteb, nstlim, ntave)
    DebugLogger.write_integration_algorithm_to_file(filename, integrator)

    [group, group, integrator] = create_upper_non_bonded_boost_integrator(system, temperature, dt, ntcmdprep, ntcmd,
                                                                          ntebprep,
                                                                          nteb, nstlim, ntave)
    DebugLogger.write_integration_algorithm_to_file(filename, integrator)

    [group, group, integrator] = create_lower_dual_non_bonded_dihederal_boost_integrator(system, temperature, dt,
                                                                                         ntcmdprep, ntcmd,
                                                                                         ntebprep, nteb, nstlim, ntave)
    DebugLogger.write_integration_algorithm_to_file(filename, integrator)

    [group, group, integrator] = create_upper_dual_non_bonded_dihederal_boost_integrator(system, temperature, dt,
                                                                                         ntcmdprep, ntcmd,
                                                                                         ntebprep, nteb, nstlim, ntave)
    DebugLogger.write_integration_algorithm_to_file(filename, integrator)

    sys.exit(22)


def create_output_directories(directories):
    for dir in directories:
        os.makedirs(dir, 0o755)


def usage():
    print("run-test.py boost-type [output-directory] [platform | device] [device] [debug] [quick] [print-algorithms]\n")
    print("\tboost-type:\t\tgamd-cmd-base|lower-total|upper-total|lower-dihedral|upper-dihedral|\n"
           "\t\t\t\tlower-dual|upper-dual|lower-nonbonded|upper-nonbonded|\n"
           "\t\t\t\tlower-dual-nonbonded-dihedral|upper-dual-nonbonded-dihedral\n")
    print("\toutput-directory:\tDirectory to output files. [default: output]\n")
    print("\tplatform:\t\tCUDA|OpenCL|CPU [defaults to OpenMM best guess for fastest]\n")
    print("\tdevice:\t\t\tUsed to specify the device index, when multiple GPUs exist on \n\t\t\t\tthe system. [defaults to CUDA, if platform not specified.]\n")
    print("\tdebug:\t\t\tpassing this string turns on debugging, which will print out system\n\t\t\t\tvariables and the integration algorithm.\n")
    print("\tquick:\t\t\tpassing this string will cause the program to run a short 60k step version\n")
    print("\tprint-algorithms:\tIntegration algorithms for all integrators will be\n\t\t\t\tsaved to all-integration-algorithms.txt in the current directory.\n")


def handle_arguments():

    quick = False
    debug = False

    if "quick" in sys.argv:
        quick = True
        sys.argv.remove("quick")

    if "debug" in sys.argv:
        debug = True
        sys.argv.remove("debug")

    if "print-algorithms" in sys.argv:
        print_integration_algorithms("all-integration-algorithms.txt")

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

    return [boost_type, output_directory, device, platform, debug, quick]


def run_simulation(unitless_temperature, dt, ntcmdprep, ntcmd, ntebprep, nteb, nstlim, ntave, boost_type,
                   output_directory, platform_name, device, running_rates: RunningRates, reweighting_offset=0,
                   debug=False):
    #pdb_filename = './data/rna_ds_separated25A_wet.pdb'
    pdb_filename = './data/rna_ds_separated25A_dry.pdb'
    starttime = time.time()
    restarting = False
    restart_checkpoint_filename = os.path.join(output_directory, "gamd.backup")
    temperature = unitless_temperature * kelvin
    pdb = PDBFile(pdb_filename)
    pdb_parmed_for_box = parmed.load_file(pdb_filename)
    #forcefield = ForceField('amoeba2018.xml')
    forcefield = ForceField('amoeba2018.xml', 'amoeba2018_gk.xml') # for implicit solvent
    #system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=HBonds)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    
    if boost_type == "upper-total" or boost_type == "upper-nonbonded":
        integrator_information = GamdIntegratorFactory.get_integrator(boost_type, system, temperature, dt, ntcmdprep,
                                                                      ntcmd, ntebprep, nteb, nstlim, ntave,
                                                                      2.5 * kilocalories_per_mole)
    elif boost_type == "upper-dual":
        integrator_information = GamdIntegratorFactory.get_integrator(boost_type, system, temperature, dt, ntcmdprep,
                                                                      ntcmd, ntebprep, nteb, nstlim, ntave,
                                                                      2.0 * kilocalories_per_mole,
                                                                      2.0 * kilocalories_per_mole)
    else:
        try:
            integrator_information = GamdIntegratorFactory.get_integrator(boost_type, system, temperature, dt,
                                                                          ntcmdprep,
                                                                          ntcmd, ntebprep, nteb, nstlim, ntave)
        except ValueError as e:
            print(e)
            usage()
            sys.exit(1)

    first_boost_group = integrator_information[0]
    second_boost_group = integrator_information[1]
    integrator = integrator_information[2]
    first_boost_type = integrator_information[3]
    second_boost_type = integrator_information[4]

    if not restarting:
        create_output_directories([output_directory, output_directory + "/states/", output_directory + "/positions/",
                                   output_directory + "/pdb/", output_directory + "/checkpoints"])

    properties = {}
    if platform_name == "CUDA":
        platform = Platform.getPlatformByName(platform_name)
        properties['CudaPrecision'] = 'mixed'
        properties['DeviceIndex'] = device
        simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    elif platform_name == "OpenCL":
        platform = Platform.getPlatformByName(platform_name)
        properties['DeviceIndex'] = device
        simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    else:
        simulation = Simulation(pdb.topology, system, integrator)

    """ # uncomment to view integrator computations
    for i in range(integrator.getNumComputations()):
        print(integrator.getComputationStep(i))

    print("exiting...")
    exit()
    """

    if restarting:
        simulation.loadCheckpoint(restart_checkpoint_filename)
        write_mode = "a"
        start_step = running_rates.get_restart_step(integrator)
        print("restarting from saved checkpoint:", restart_checkpoint_filename,
              "at step:", start_step)
        running_range = running_rates.get_restart_batch_run_range(integrator)
    else:
        running_range = running_rates.get_batch_run_range()
        
        """
        simulation.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        """
        simulation.context.setPositions(pdb.positions)
        box_vectors = pdb_parmed_for_box.box_vectors
        simulation.context.setPeriodicBoxVectors(*box_vectors)
        simulation.minimizeEnergy()
        #
        # Should we be setting this to the temperature?
        #
        simulation.context.setVelocitiesToTemperature(unitless_temperature * kelvin)
        simulation.saveState(output_directory + "/states/initial-state.xml")
        simulation.reporters.append(DCDReporter(output_directory + '/output.dcd', running_rates.get_save_rate()))
        simulation.reporters.append(
            utils.ExpandedStateDataReporter(system,
                                            output_directory + '/state-data.log',
                                            running_rates.get_reporting_rate(),
                                            step=True,
                                            brokenOutForceEnergies=True, temperature=True, potentialEnergy=True,
                                            totalEnergy=True, volume=True))
        print("startup time: \t", time.time() - starttime)
        write_mode = "w"
    gamd_log_filename = os.path.join(output_directory, "gamd.log")
    gamd_reweighting_filename = os.path.join(output_directory, "gamd-reweighting.log")

    gamd_logger = GamdLogger(gamd_log_filename, write_mode, integrator, simulation,
                             first_boost_type, first_boost_group,
                             second_boost_type, second_boost_group)

    gamd_reweighting_logger = GamdLogger(gamd_reweighting_filename, write_mode, integrator, simulation,
                                         first_boost_type, first_boost_group,
                                         second_boost_type, second_boost_group)

    start_production_logging_step = ntcmd + nteb + (running_rates.get_batch_run_rate() * reweighting_offset)

    with open(output_directory + "/" + "production-start-step.txt", "w") as prodstartstep_file:
        prodstartstep_file.write(str(start_production_logging_step))

    if not restarting:
        gamd_logger.write_header()
        gamd_reweighting_logger.write_header()

    print("Running: \t " + str(integrator.get_total_simulation_steps()) + " steps")

    if debug:
        ignore_fields = {"stageOneIfValueIsZeroOrNegative", "stageTwoIfValueIsZeroOrNegative",
                         "stageThreeIfValueIsZeroOrNegative", "stageFourIfValueIsZeroOrNegative",
                         "stageFiveIfValueIsZeroOrNegative", "thermal_energy", "collision_rate",
                         "vscale", "fscale", "noisescale"}
        debug_filename = os.path.join(output_directory, "debug.csv")
        debug_logger = DebugLogger(debug_filename, write_mode, ignore_fields)
        print("Debugging enabled.")
        debug_logger.write_integration_algorithm_to_file(output_directory + "/integration-algorithm.txt", integrator)
        debug_logger.write_global_variables_headers(integrator)

    start_date_time = datetime.datetime.now()
    batch_run_rate = running_rates.get_batch_run_rate()
    for batch_frame in running_range:
        step = running_rates.get_step_from_frame(batch_frame)

        if running_rates.is_save_step(step):
            simulation.saveCheckpoint(restart_checkpoint_filename)

        # TEST
        #            if step == 250 and not restarting:
        #                print("sudden, unexpected interruption!")
        #                exit()
        # END TEST

        gamd_logger.mark_energies()
        if running_rates.is_save_step(step):
            gamd_reweighting_logger.mark_energies()

        try:

            #
            #  NOTE:  We need to save off the starting total and dihedral potential energies, since we
            #         end up modifying them by updating the state of the particles.  This allows us to write
            #         out the values as they were at the beginning of the step for what all of the calculations
            #         for boosting were based on.
            #

            simulation.step(batch_run_rate)
            if debug and running_rates.is_debugging_step(batch_frame):
                debug_logger.write_global_variables_values(integrator)

            if running_rates.is_save_step(step):
                gamd_logger.write_to_gamd_log(step)
                if step >= start_production_logging_step:
                    gamd_reweighting_logger.write_to_gamd_log(step)

        except Exception as e:
            print("Failure on step " + str(step))
            print(e)
            if debug:
                debug_logger.print_global_variables_to_screen(integrator)
                debug_logger.write_global_variables_values(integrator)
            gamd_logger.write_to_gamd_log(step)
            gamd_reweighting_logger.write_to_gamd_log(step)

            sys.exit(2)

        if (step % ntave) == 0:
            simulation.saveState(output_directory + "/states/" + str(step) + ".xml")
            step_checkpoint_filename = os.path.join(output_directory, "checkpoints", str(step) + ".bin")
            positions_filename = os.path.join(output_directory, "positions",
                                                     "coordinates-" + str(step) + ".csv")
            simulation.saveCheckpoint(step_checkpoint_filename)
            integrator.create_positions_file(positions_filename)

    return start_date_time


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
