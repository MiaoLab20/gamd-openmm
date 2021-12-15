"""
Run a GaMD calculation

More here...
"""

import os
import sys
import time
import traceback
import pprint
import argparse
import glob
import datetime
import shutil
import subprocess

# TODO: need to remove import * - does not conform to PEP8
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from gamd import utils as utils
from gamd import parser
from gamd import gamdSimulation
from gamd.stage_integrator import BoostType
from gamd.DebugLogger import DebugLogger
from gamd.GamdLogger import GamdLogger

RESTART_CHECKPOINT_FILENAME = "gamd_restart.checkpoint"

def create_output_directories(directories, overwrite_output=False):
    if overwrite_output:
        print("Overwrite output set to True")
        for dir in directories:
            if os.path.exists(dir):
                print("Deleting old output directory:", dir)
                os.system('rm -r %s' % dir)
    
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

class RunningRates:

    def __init__(self, number_of_simulation_steps: int, save_rate: int, reporting_rate: int,
                 debugging_enabled: bool=False, debugging_step_function=None) -> None:
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

class Runner:
    def __init__(self, config, gamdSimulation, debug):
        self.config = config
        self.gamdSim = gamdSimulation
        self.debug = debug
        nstlim = self.config.integrator.number_of_steps.total_simulation_length
        self.chunk_size = self.config.outputs.reporting.compute_chunk_size()
        if debug:
            self.running_rates = RunningRates(nstlim, self.chunk_size, 1, True)
        else:
            self.running_rates = RunningRates(nstlim, self.chunk_size, 
                                              self.chunk_size, False)
        return
    
    def run(self, restart=False):
        debug = self.debug
        start_date_time = datetime.datetime.now()
        print("Start Time: \t", start_date_time.strftime("%b-%d-%Y    %H:%M"))
        initial_temperature = self.config.temperature
        target_temperature = self.config.temperature
        output_directory = self.config.outputs.directory
        restart_checkpoint_interval \
            = self.config.outputs.reporting.restart_checkpoint_interval
        statistics_interval = self.config.outputs.reporting.statistics_interval
        restart_checkpoint_filename = os.path.join(
            output_directory, RESTART_CHECKPOINT_FILENAME)
        chunk_size = self.chunk_size
        overwrite_output = self.config.outputs.overwrite_output
        if not restart:
            create_output_directories(
                [output_directory, 
                 os.path.join(output_directory, "states/"), 
                 os.path.join(output_directory, "positions/"),
                 os.path.join(output_directory, "pdb/"),
                 os.path.join(output_directory, "checkpoints")],
                overwrite_output)
        
        system = self.gamdSim.system
        simulation = self.gamdSim.simulation
        integrator = self.gamdSim.integrator
        traj_reporter = self.gamdSim.traj_reporter
        ntcmdprep = self.config.integrator.number_of_steps.conventional_md_prep
        ntcmd = self.config.integrator.number_of_steps.conventional_md
        ntebprep = self.config.integrator.number_of_steps.gamd_equilibration_prep
        nteb = self.config.integrator.number_of_steps.gamd_equilibration
        nstlim = self.config.integrator.number_of_steps.total_simulation_length
        ntave = self.config.integrator.number_of_steps.averaging_window_interval
        #group = self.config.dihedral_group
        
        extension = self.config.outputs.reporting.coordinates_file_type
        if restart:
            simulation.loadCheckpoint(restart_checkpoint_filename)
            currentStep = int(integrator.getGlobalVariableByName("stepCount"))
            simulation.currentStep = currentStep
            write_mode = "a"
            start_chunk = (currentStep // chunk_size) + 1
            print("restarting from saved checkpoint:", 
                  restart_checkpoint_filename, "at step:", start_chunk)
            # see how many restart files have already been created
            state_data_restart_files_glob = os.path.join(
                output_directory, 'state-data.restart*.log')
            state_data_restarts_list = glob.glob(state_data_restart_files_glob)
            restart_index = len(state_data_restarts_list) + 1
            state_data_name = os.path.join(
                output_directory, 'state-data.restart%d.log' % restart_index)
            traj_name = os.path.join(
                self.config.outputs.directory, 
                'output.restart%d.%s' % (restart_index, extension))
            running_range = self.running_rates.get_restart_batch_run_range(integrator)
            
        else:
            running_range = self.running_rates.get_batch_run_range()
            simulation.saveState(
                os.path.join(output_directory, "states/initial-state.xml"))
            write_mode = "w"
            start_chunk = 1
            traj_name = os.path.join(
                self.config.outputs.directory, 'output.%s' % extension)
            state_data_name = os.path.join(output_directory, 'state-data.log')
            
        if traj_reporter:
            simulation.reporters.append(traj_reporter(
                traj_name, self.config.outputs.reporting.coordinates_interval))
        simulation.reporters.append(utils.ExpandedStateDataReporter(
            system, state_data_name, 
            self.config.outputs.reporting.energy_interval, step=True, 
            brokenOutForceEnergies=True, temperature=True, 
            potentialEnergy=True, totalEnergy=True, 
            volume=True))
            
        
        # TODO: check if we really want to get this quantity from integrator
        # instead of the config object
        end_chunk = int(integrator.get_total_simulation_steps() \
                        // chunk_size) + 1
        
        #####################################
        # NEW CODE STARTED HERE
        #####################################
        reweighting_offset = 0
        gamd_log_filename = os.path.join(output_directory, "gamd.log")
        gamd_reweighting_filename = os.path.join(output_directory, "gamd-reweighting.log")
        gamd_logger = GamdLogger(gamd_log_filename, write_mode, integrator, simulation,
                                 self.gamdSim.first_boost_type, self.gamdSim.first_boost_group,
                                 self.gamdSim.second_boost_type, self.gamdSim.second_boost_group)
    
        gamd_reweighting_logger = GamdLogger(gamd_reweighting_filename, write_mode, integrator, simulation,
                                             self.gamdSim.first_boost_type, self.gamdSim.first_boost_group,
                                             self.gamdSim.second_boost_type, self.gamdSim.second_boost_group)
        
        start_production_logging_step = ntcmd + nteb + (self.running_rates.get_batch_run_rate() * reweighting_offset)
    
        with open(output_directory + "/" + "production-start-step.txt", "w") as prodstartstep_file:
            prodstartstep_file.write(str(start_production_logging_step))
    
        if not restart:
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
        batch_run_rate = self.running_rates.get_batch_run_rate()
        for batch_frame in running_range:
            step = self.running_rates.get_step_from_frame(batch_frame)
    
            if self.running_rates.is_save_step(step):
                simulation.saveCheckpoint(restart_checkpoint_filename)
    
            # TEST
            #            if step == 250 and not restarting:
            #                print("sudden, unexpected interruption!")
            #                exit()
            # END TEST
    
            gamd_logger.mark_energies()
            if self.running_rates.is_save_step(step):
                gamd_reweighting_logger.mark_energies()
    
            try:
    
                #
                #  NOTE:  We need to save off the starting total and dihedral potential energies, since we
                #         end up modifying them by updating the state of the particles.  This allows us to write
                #         out the values as they were at the beginning of the step for what all of the calculations
                #         for boosting were based on.
                #
    
                simulation.step(batch_run_rate)
                if debug and self.running_rates.is_debugging_step(batch_frame):
                    debug_logger.write_global_variables_values(integrator)
    
                if self.running_rates.is_save_step(step):
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
            
        """
        with open(os.path.join(output_directory, "gamd.log"), write_mode) \
                as gamdLog:
            if not restart:
                gamdLog.write(
                    "# Gaussian accelerated Molecular Dynamics log file\n")
                gamdLog.write(
                    "# All energy terms are stored in unit of kcal/mol\n")
                gamdLog.write(
                    "# ntwx,total_nstep,Unboosted-Potential-Energy,"\
                    "Unboosted-Dihedral-Energy,Total-Force-Weight,"\
                    "Dihedral-Force-Weight,Boost-Energy-Potential,"\
                    "Boost-Energy-Dihedral\n")
            
            for chunk in range(start_chunk, end_chunk):
                if chunk % (restart_checkpoint_interval // chunk_size) == 0:
                    simulation.saveCheckpoint(restart_checkpoint_filename)
                
                '''
                # TEST
                if chunk == 249:
                    print("sudden, unexpected interruption!")
                    exit()
                # END TEST
                '''
                if chunk % (statistics_interval // chunk_size) == 0:
                    state = simulation.context.getState(getEnergy=True)
                    dihedral_state = simulation.context.getState(getEnergy=True, groups={group})
                    force_scaling_factors = integrator.get_force_scaling_factors()
                    boost_potentials = integrator.get_boost_potentials()
                    total_potential_energy = str(state.getPotentialEnergy() / (kilojoules_per_mole * 4.184))
                    dihedral_energy = str(dihedral_state.getPotentialEnergy() / (kilojoules_per_mole * 4.184))
                    total_force_scaling_factor = str(force_scaling_factors["ForceScalingFactor_" + BoostType.TOTAL.value])
                    dihedral_force_scaling_factor = str(force_scaling_factors["ForceScalingFactor_Dihedral"])
                    total_boost_potential = str(boost_potentials["BoostPotential_" + BoostType.TOTAL.value] / 4.184)
                    dihedral_boost_potential = str(boost_potentials["BoostPotential_Dihedral"] / 4.184)
                
                try:
                    simulation.step(chunk_size)
                    #print_global_variables(integrator)
                    state = simulation.context.getState(
                        getEnergy=True, groups={group})
                    
                    # TODO: deal with these strange units issues
                    gamdLog.write("\t".join((
                        str(chunk_size), str(chunk*chunk_size), 
                        total_potential_energy, dihedral_energy,
                        total_force_scaling_factor, 
                        dihedral_force_scaling_factor,
                        total_boost_potential, dihedral_boost_potential))+"\n")
    
                except Exception as e:
                    print("Failure on step " + str(chunk*chunk_size))
                    print(e)
                    state = simulation.context.getState(
                        getEnergy=True, groups={group})
                    gamdLog.write("\t".join((
                        "Fail Step: " + str(chunk_size), str(chunk*chunk_size), 
                        total_potential_energy, dihedral_energy,
                        total_force_scaling_factor, 
                        dihedral_force_scaling_factor,
                        total_boost_potential, dihedral_boost_potential))+"\n")
                    sys.exit(2)
                
                assert integrator.ntave >= chunk_size, "The number of steps "\
                    "per averaging must be greater than or equal to chunk_size."
                if chunk % (integrator.ntave // chunk_size) == 0:
    
                    simulation.saveState(
                        os.path.join(output_directory, 
                                     "states", str(chunk*chunk_size) + ".xml"))
                    simulation.saveCheckpoint(
                        os.path.join(output_directory, "checkpoints",
                                     str(chunk*chunk_size) + ".bin"))
                    positions_filename = os.path.join(output_directory,
                        "positions", "coordinates-" + str(chunk*chunk_size) \
                        + ".csv")
                    integrator.create_positions_file(positions_filename)
        """
        
        print("Start Time: \t", start_date_time.strftime("%b-%d-%Y    %H:%M"))
        end_date_time = datetime.datetime.now()
        print("End Time: \t", end_date_time.strftime("%b-%d-%Y    %H:%M"))
    
        time_difference = end_date_time - start_date_time
        steps_per_second = nstlim / time_difference.seconds
        print("Execution rate for this run:  ", str(steps_per_second), " steps per second.")
        daily_execution_rate = steps_per_second * 3600 * 24 * self.config.integrator.dt
        print("Daily execution rate:         ", str(daily_execution_rate.value_in_unit(nanoseconds)), " ns per day.")
        production_starting_frame = ((ntcmd + nteb) / chunk_size) + reweighting_offset
        run_post_simulation(self.config.temperature, output_directory, production_starting_frame)


def main():
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "input_file_type", metavar="INPUT_FILE_TYPE", type=str, 
        help="The type of file being provided. Available options are: 'xml', "\
        "... More to come later")
    argparser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, 
        help="name of input file for GaMD calculation. only XML format is "\
            "currently preferred.")
    argparser.add_argument("-r", "--restart", dest="restart", default=False,
                           help="Restart simulation from backup checkpoint in "\
                           "input file", action="store_true")
    argparser.add_argument("-p", "--platform", dest="platform", default="CUDA",
                           help="Define the platform that will run the "\
                           "simulations. Default is 'CUDA', but other options"\
                           "include: 'reference', 'CPU', and 'OpenCL'.", 
                           type=str)
    argparser.add_argument("-d", "--device_index", dest="device_index", default="0",
                           help="modify which device_index to run the "\
                           "simulation on. For example, the number 0 or 1 "\
                           "would suffice. To run on multiple GPU indices, "
                           "simply enter comma separated indices. Example: "\
                           "'0,1'. If a value is not supplied, the value '0' "\
                           "will be used by default.", type=str)
    argparser.add_argument("-D", "--debug", dest="debug", default=False,
                           help="Whether to start the run in debug mode.", 
                           action="store_true")
    
    args = argparser.parse_args()  # parse the args into a dictionary
    args = vars(args)
    config_file_type = args["input_file_type"]
    config_filename = args["input_file"]
    restart = args["restart"]
    platform = args["platform"]
    cuda_index = args["device_index"]
    debug = args["debug"]
    
    parserFactory = parser.ParserFactory()
    config = parserFactory.parse_file(config_filename, config_file_type)
    gamdSimulationFactory = gamdSimulation.GamdSimulationFactory()
    gamdSim = gamdSimulationFactory.createGamdSimulation(
        config, platform, cuda_index)
    
    # If desired, modify OpenMM objects in gamdSimulation object here...
    
    runner = Runner(config, gamdSim, debug)
    runner.run(restart)
    
               
if __name__ == "__main__":
    main()
