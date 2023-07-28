import time 
import datetime
import glob
import os
import shutil
import subprocess
import sys

import mdtraj
import openmm.unit as unit
import openmm.app as openmm_app

from gamd import utils as utils
from gamd.DebugLogger import DebugLogger, NoOpDebugLogger
from gamd.GamdLogger import GamdLogger, NoOpGamdLogger


def create_output_directories(directories, overwrite_output=False):
    if overwrite_output:
        for directory in directories:
            if os.path.exists(directory):
                print("Deleting old output directory:", directory)
                os.system('rm -r %s' % directory)
    
    for directory in directories:
        os.makedirs(directory, 0o755)


def get_global_variable_names(integrator):
    for index in range(0, integrator.getNumGlobalVariables()):
        print(integrator.getGlobalVariableName(index))


def print_global_variables(integrator):
    for index in range(0, integrator.getNumGlobalVariables()):
        name = integrator.getGlobalVariableName(index)
        value = integrator.getGlobalVariableByName(name)
        print(name + ":  " + str(value))


def write_gamd_production_restart_file(output_directory, integrator,
                                       first_boost_type, second_boost_type):
    gamd_prod_restart_filename = os.path.join(output_directory, "gamd-restart.dat")
    values = integrator.get_statistics()

    with open(gamd_prod_restart_filename, "w") as gamd_prod_restart_file:
        for key in values.keys():
            gamd_prod_restart_file.write(key + "=" + str(values[key]) + "\n")


def get_config_and_simulation_values(gamd_simulation, config):
    output_directory = config.outputs.directory
    overwrite_output = config.outputs.overwrite_output
    system = gamd_simulation.system
    simulation = gamd_simulation.simulation
    integrator = gamd_simulation.integrator
    dt = config.integrator.dt
    ntcmdprep = config.integrator.number_of_steps.conventional_md_prep
    ntcmd = config.integrator.number_of_steps.conventional_md
    ntebprep = config.integrator.number_of_steps.gamd_equilibration_prep
    nteb = config.integrator.number_of_steps.gamd_equilibration
    last_step_of_equilibration = ntcmd + nteb
    nstlim = config.integrator.number_of_steps.total_simulation_length
    ntave = config.integrator.number_of_steps.averaging_window_interval

    return [output_directory, overwrite_output, system, simulation, dt,
            integrator, ntcmdprep, ntcmd, ntebprep,
            nteb, last_step_of_equilibration, nstlim, ntave]


def print_runtime_information(start_date_time, dt, nstlim, current_step):
    end_date_time = datetime.datetime.now()
    time_difference = end_date_time - start_date_time
    if time_difference.seconds > 0:
        steps_per_second = (nstlim - current_step) / time_difference.seconds
    else:
        steps_per_second = 0
    daily_execution_rate = (steps_per_second * 3600 * 24 * dt)

    print("Start Time: \t", start_date_time.strftime("%b-%d-%Y    %H:%M:%S"))
    print("End Time: \t", end_date_time.strftime("%b-%d-%Y    %H:%M:%S"))
    print("Execution rate for this run:  ", str(steps_per_second),
          " steps per second.")
    print("Daily execution rate:         ",
          str(daily_execution_rate.value_in_unit(unit.nanoseconds)),
          " ns per day.")


class RunningRates:

    def __init__(self, number_of_simulation_steps: int, save_rate: int,
                 reporting_rate: int,
                 debugging_enabled: bool=False,
                 debugging_step_function=None) -> None:
        """Constructor

            The save_rate determines that rate at which to write out
            DCD/PDB file, checkpoint, coordinates, and write to the GaMD log
            files. The value of save_rate should evenly divide the number of
            simulation steps.

            The reporting_rate determines the rate at which to write out to
            the state-data.log and the debugging files, if enabled.  The value
            of the reporting rate should evenly divide the number of simulation
             steps.

            The save_rate or reporting_rate should be a multiple of the other
            value, if debugging is enabled.

            The save_rate, reporting_rate, and debugging flag determine the
            batch_run_rate, which is just the number of steps to have
            OpenMM execute at one time in the simulation.

        Parameters

        :param number_of_simulation_steps:  (int) Same as nstlim, the total
                                            number of steps in the simulation.

        :param save_rate: (int) the number of steps before performing a save.

        :param reporting_rate:  (int) the number of steps before printing to
                                reports.

        :param debugging_enabled:   (boolean) indicates whether debugging is
                                    enabled.  Defaults: False

        """
        if ((number_of_simulation_steps % save_rate) != 0 and
                (number_of_simulation_steps % reporting_rate) != 0):
            raise ValueError("RunningRates:  The save_rate and "
                             "reporting_rate should evenly divide "
                             "into the number of steps in the simulation.")

        if debugging_enabled and not (((save_rate % reporting_rate) == 0) or
                                      ((reporting_rate % save_rate) == 0)):
            raise ValueError("RunningRates:  When debugging is enabled, "
                             "save_rate/reporting_rate must be a multiple of "
                             "the other value.")

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
        return range(1, self.number_of_simulation_steps //
                     self.get_batch_run_rate()+1)

    def get_restart_step(self, integrator):
        start_step = int(round(integrator.getGlobalVariableByName("stepCount") /
                         self.get_batch_run_rate()))
        return start_step

    def get_restart_batch_run_range(self, integrator):
        return range(self.get_restart_step(integrator)+1,
                     self.number_of_simulation_steps //
                     self.get_batch_run_rate()+1)

    def get_step_from_frame(self, frame):
        return frame * self.get_batch_run_rate()


class Runner:
    def __init__(self, config, gamd_simulation, debug):
        self.config = config
        self.gamd_simulation = gamd_simulation
        self.debug = debug
        nstlim = self.config.integrator.number_of_steps.total_simulation_length
        self.chunk_size = self.config.outputs.reporting.compute_chunk_size()
        if debug:
            self.running_rates = RunningRates(nstlim, self.chunk_size, 1, True)
        else:
            self.running_rates = RunningRates(nstlim, self.chunk_size, 
                                              self.chunk_size, False)
        self.gamd_logger_enabled = True
        self.gamd_reweighting_logger_enabled = False
        self.state_data_reporter_enabled = False
        self.gamd_dat_reporter_enabled = False
        return

    def run_post_simulation(self, temperature, output_directory,
                            production_starting_frame):
        return

    def save_initial_configuration(self, production_logging_start_step,
                                   temperature):

        config_filename = os.path.join(self.config.outputs.directory,
                                       "input.xml")
        prodstartstep_filename = os.path.join(self.config.outputs.directory,
                                              "production-start-step.txt")
        temperature_filename = os.path.join(self.config.outputs.directory,
                                            "temperature.dat")

        self.config.serialize(config_filename)

        with open(temperature_filename, "w") as temperature_file:
            temperature_file.write(str(temperature))

        with open(prodstartstep_filename, "w") as prodstartstep_file:
            prodstartstep_file.write(str(production_logging_start_step))

    def register_trajectory_reporter(self, restart):
        simulation = self.gamd_simulation.simulation
        traj_reporter = self.gamd_simulation.traj_reporter
        output_directory = self.config.outputs.directory
        extension = self.config.outputs.reporting.coordinates_file_type
        traj_name = os.path.join(output_directory, 'output.%s' % extension)
        traj_append = restart

        if traj_reporter == openmm_app.DCDReporter:
            simulation.reporters.append(traj_reporter(
                traj_name, self.config.outputs.reporting.coordinates_interval,
                append=traj_append))
        elif traj_reporter == openmm_app.PDBReporter:
            simulation.reporters.append(traj_reporter(
                traj_name, self.config.outputs.reporting.coordinates_interval))
        elif traj_reporter == mdtraj.reporters.HDF5Reporter:
            simulation.reporters.append(traj_reporter(
                traj_name, self.config.outputs.reporting.coordinates_interval))

    def register_state_data_reporter(self, restart):
        if self.state_data_reporter_enabled:
            simulation = self.gamd_simulation.simulation
            output_directory = self.config.outputs.directory
            system = self.gamd_simulation.system

            if restart:
                state_data_restart_files_glob = os.path.join(
                    output_directory, 'state-data.restart*.log')
                state_data_restarts_list = glob.glob(state_data_restart_files_glob)
                restart_index = len(state_data_restarts_list) + 1
                state_data_name = os.path.join(
                    output_directory, 'state-data.restart%d.log' % restart_index)
            else:
                state_data_name = os.path.join(output_directory, 'state-data.log')

            simulation.reporters.append(utils.ExpandedStateDataReporter(
                system, state_data_name,
                self.config.outputs.reporting.energy_interval, step=True,
                brokenOutForceEnergies=True, temperature=True,
                potentialEnergy=True, totalEnergy=True,
                volume=True))

    def register_gamd_data_reporter(self, restart):
        if self.gamd_dat_reporter_enabled:
            simulation = self.gamd_simulation.simulation
            output_directory = self.config.outputs.directory
            gamd_running_dat_filename = os.path.join(output_directory, "gamd-running.csv")
            integrator = self.gamd_simulation.integrator

            #
            #  The GamdDatReporter uses the write mode to determine whether to write out headers or not.
            #

            if restart:
                write_mode = "a"
            else:
                write_mode = "w"
            gamd_dat_reporter = utils.GamdDatReporter(gamd_running_dat_filename, write_mode,
                                                      integrator)
            simulation.reporters.append(gamd_dat_reporter)

    def register_gamd_logger(self, restart):
        if self.gamd_logger_enabled:
            output_directory = self.config.outputs.directory
            gamd_log_filename = os.path.join(output_directory, "gamd.log")
            integrator = self.gamd_simulation.integrator
            simulation = self.gamd_simulation.simulation
            if restart:
                write_mode = "a"
            else:
                write_mode = "w"

            gamd_logger = GamdLogger(gamd_log_filename, write_mode, integrator,
                                     simulation, self.gamd_simulation.first_boost_type,
                                     self.gamd_simulation.first_boost_group,
                                     self.gamd_simulation.second_boost_type,
                                     self.gamd_simulation.second_boost_group)
            if not restart:
                gamd_logger.write_header()
        else:
            gamd_logger = NoOpGamdLogger()
        return gamd_logger

    def register_gamd_reweighting_logger(self, restart):
        if self.gamd_reweighting_logger_enabled:
            output_directory = self.config.outputs.directory
            integrator = self.gamd_simulation.integrator
            simulation = self.gamd_simulation.simulation
            gamd_reweighting_filename = os.path.join(output_directory,
                                                     "gamd-reweighting.log")
            if restart:
                write_mode = "a"
            else:
                write_mode = "w"

            gamd_reweighting_logger = GamdLogger(gamd_reweighting_filename,
                                                 write_mode, integrator,
                                                 simulation,
                                                 self.gamd_simulation.first_boost_type,
                                                 self.gamd_simulation.first_boost_group,
                                                 self.gamd_simulation.second_boost_type,
                                                 self.gamd_simulation.second_boost_group)
            if not restart:
                gamd_reweighting_logger.write_header()
        else:
            gamd_reweighting_logger = NoOpGamdLogger()
        return gamd_reweighting_logger

    def register_debug_logger(self, restart):
        output_directory = self.config.outputs.directory
        integrator = self.gamd_simulation.integrator
        debug_filename = os.path.join(output_directory, "debug.csv")
        if restart:
            write_mode = "a"
        else:
            write_mode = "w"

        if self.debug:
            ignore_fields = {"stageOneIfValueIsZeroOrNegative",
                             "stageTwoIfValueIsZeroOrNegative",
                             "stageThreeIfValueIsZeroOrNegative",
                             "stageFourIfValueIsZeroOrNegative",
                             "stageFiveIfValueIsZeroOrNegative",
                             "thermal_energy", "collision_rate",
                             "vscale", "fscale", "noisescale"}

            debug_logger = DebugLogger(debug_filename, write_mode,
                                       ignore_fields)
            print("Debugging enabled.")
            int_algorithm_filename = os.path.join(output_directory,
                                                  "integration-algorithm.txt")
            debug_logger.write_integration_algorithm_to_file(
                int_algorithm_filename, integrator)

            debug_logger.write_global_variables_headers(integrator)
        else:
            debug_logger = NoOpDebugLogger()

        return debug_logger

    def run(self, restart=False):
        chunk_size = self.chunk_size
        output_directory, overwrite_output, system, simulation, dt, \
            integrator, ntcmdprep, ntcmd, ntebprep, nteb, \
            last_step_of_equilibration, nstlim, ntave \
            = get_config_and_simulation_values(self.gamd_simulation, self.config)

        restart_checkpoint_filename = os.path.join(
            output_directory, "gamd_restart.checkpoint")

        if not restart:
            create_output_directories(
                [output_directory],
                overwrite_output)

        if restart:
            simulation.loadCheckpoint(restart_checkpoint_filename)
            state = simulation.context.getState()
            state_time = state.getTime().value_in_unit(unit.picoseconds)
            integrator_dt = dt.value_in_unit(unit.picoseconds)
            current_step = int(round(state_time / integrator_dt))
            simulation.currentStep = current_step
            print("restarting from saved checkpoint:",
                  restart_checkpoint_filename, "at step:", current_step)

            running_range = self.running_rates.get_restart_batch_run_range(
                integrator)
        else:
            current_step = 0
            running_range = self.running_rates.get_batch_run_range()

        self.register_trajectory_reporter(restart)
        self.register_state_data_reporter(restart)
        self.register_gamd_data_reporter(restart)
        debug_logger = self.register_debug_logger(restart)
        gamd_logger = self.register_gamd_logger(restart)
        gamd_reweighting_logger = self.register_gamd_reweighting_logger(restart)

        reweighting_offset = 0
        production_logging_start_step = (ntcmd + nteb +
                                         (self.running_rates.
                                          get_batch_run_rate()
                                          * reweighting_offset))

        self.save_initial_configuration(production_logging_start_step,
                                        self.config.temperature)

        print("Running: \t ",
              str(integrator.get_total_simulation_steps() - current_step),
              " steps")

        start_date_time = datetime.datetime.now()
        start_time = time.time()
        batch_run_rate = self.running_rates.get_batch_run_rate()
        for batch_frame in running_range:
            step = self.running_rates.get_step_from_frame(batch_frame)

            if self.running_rates.is_save_step(step):
                simulation.saveCheckpoint(restart_checkpoint_filename)

            gamd_logger.mark_energies()
            if self.running_rates.is_save_step(step):
                gamd_reweighting_logger.mark_energies()

            try:

                #
                #  NOTE:  We need to save off the starting total and dihedral
                #  potential energies, since we end up modifying them by
                #  updating the state of the particles.  This allows us to
                #  write out the values as they were at the beginning of the
                #  step for what all of the calculations for boosting were
                #  based on.
                #

                simulation.step(batch_run_rate)
                if self.running_rates.is_debugging_step(batch_frame):
                    debug_logger.write_global_variables_values(integrator)

                if self.running_rates.is_save_step(step):
                    gamd_logger.write_to_gamd_log(step)
                    if step >= production_logging_start_step:
                        gamd_reweighting_logger.write_to_gamd_log(step)

            except Exception as e:
                print("Failure on step " + str(step))
                print(e)
                gamd_logger.close()
                gamd_reweighting_logger.close()
                debug_logger.print_global_variables_to_screen(integrator)
                debug_logger.write_global_variables_values(integrator)
                debug_logger.close()

                sys.exit(2)

            if step == last_step_of_equilibration:
                write_gamd_production_restart_file(output_directory, integrator,
                                                   self.gamd_simulation.first_boost_type,
                                                   self.gamd_simulation.second_boost_type)

        #
        # These calls are here to guarantee that the file buffers have been
        # flushed, prior to any post-simulations steps attempting
        # to utilize these files.
        #
        gamd_logger.close()
        gamd_reweighting_logger.close()
        debug_logger.close()

        simulation.saveCheckpoint(restart_checkpoint_filename)
        print_runtime_information(start_date_time, dt, nstlim, current_step)
        production_starting_frame = (((ntcmd + nteb) / chunk_size) +
                                     reweighting_offset)
        self.run_post_simulation(self.config.temperature, output_directory,
                                 production_starting_frame)


class DeveloperRunner(Runner):
    def __init__(self, config, gamd_simulation, debug):
        super().__init__(config, gamd_simulation, debug)
        self.gamd_logger_enabled = True
        self.gamd_reweighting_logger_enabled = True
        self.state_data_reporter_enabled = True
        self.gamd_dat_reporter_enabled = True
        return


class NoLogRunner(Runner):
    def __init__(self, config, gamd_simulation, debug):
        super().__init__(config, gamd_simulation, debug)
        self.gamd_logger_enabled = False
        self.gamd_reweighting_logger_enabled = False
        self.state_data_reporter_enabled = False
        self.gamd_dat_reporter_enabled = False
        return


