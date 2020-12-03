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

# TODO: need to remove import * - does not conform to PEP8
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from gamd import utils as utils
from gamd import parser
from gamd import gamdSimulation
from gamd.stage_integrator import BoostType


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

class Runner:
    def __init__(self, config, gamdSimulation):
        self.config = config
        self.gamdSim = gamdSimulation
    
    def run(self, restart=False):
        initial_temperature = self.config.initial_temperature
        target_temperature = self.config.target_temperature
        output_directory = self.config.output_directory
        restart_checkpoint_frequency = self.config.restart_checkpoint_frequency
        restart_checkpoint_filename = os.path.join(
            output_directory, self.config.restart_checkpoint_filename)
        chunk_size = self.config.chunk_size
        overwrite_output = self.config.overwrite_output
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
        #force_groups = self.gamdSim.force_groups
        traj_reporter = self.gamdSim.traj_reporter
        group = self.config.dihedral_group
        
        extension = self.config.coordinates_reporter_file_type
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
                self.config.output_directory, 
                'output.restart%d.%s' % (restart_index, extension))
            
        else:
            simulation.saveState(
                os.path.join(output_directory, "states/initial-state.xml"))
            write_mode = "w"
            start_chunk = 1
            traj_name = os.path.join(
                self.config.output_directory, 'output.%s' % extension)
            state_data_name = os.path.join(output_directory, 'state-data.log')
            
        if traj_reporter:
            simulation.reporters.append(traj_reporter(
                traj_name, self.config.coordinates_reporter_frequency,))
        simulation.reporters.append(utils.ExpandedStateDataReporter(
            system, state_data_name, 
            self.config.energy_reporter_frequency, step=True, 
            brokenOutForceEnergies=True, temperature=True, 
            potentialEnergy=True, totalEnergy=True, 
            volume=True))
            
        
        # TODO: check if we really want to get this quantity from integrator
        # instead of the config object
        end_chunk = int(integrator.get_total_simulation_steps() \
                        // chunk_size) + 1
        
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
                if chunk % (restart_checkpoint_frequency // chunk_size) == 0:
                    simulation.saveCheckpoint(restart_checkpoint_filename)
                
                '''
                # TEST
                if chunk == 249:
                    print("sudden, unexpected interruption!")
                    exit()
                # END TEST
                '''
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
                    
def main():
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        'input_file_type', metavar='INPUT_FILE_TYPE', type=str, 
        help="The type of file being provided. Available options are: 'xml', "\
        "... More to come later")
    argparser.add_argument(
        'input_file', metavar='INPUT_FILE', type=str, 
        help="name of input file for GaMD calculation. A variety of input"\
        "formats are allowed, but XML format is preferred")
    argparser.add_argument('-r', '--restart', dest='restart', default=False,
                           help="Restart simulation from backup checkpoint in "\
                           "input file", action="store_true")
    
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    config_file_type = args['input_file_type']
    config_filename = args['input_file']
    restart = args['restart']
    
    parserFactory = parser.ParserFactory()
    config = parserFactory.parse_file(config_filename, config_file_type)
    gamdSimulationFactory = gamdSimulation.GamdSimulationFactory()
    gamdSim = gamdSimulationFactory.createGamdSimulation(config)
    
    # If desired, modify OpenMM objects in gamdSimulation object here...
    
    runner = Runner(config, gamdSim)
    runner.run(restart)
    
               
if __name__ == "__main__":
    main()
