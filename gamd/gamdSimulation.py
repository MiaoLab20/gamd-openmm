"""
Created on Oct 28, 2020

Creates all OpenMM objects from Config() object that can be used in a
GaMD simulation.

@author: lvotapka
"""
import os

import parmed
import openmm as openmm
import openmm.app as openmm_app
import openmm.unit as unit

from gamd import parser
# change to generic integrator someday
from gamd.langevin.total_boost_integrators import LowerBoundIntegrator as TotalLowerBoundIntegrator
from gamd.langevin.total_boost_integrators import UpperBoundIntegrator as TotalUpperBoundIntegrator
from gamd.langevin.dihedral_boost_integrators import LowerBoundIntegrator as DihedralLowerBoundIntegrator
from gamd.langevin.dihedral_boost_integrators import UpperBoundIntegrator as DihedralUpperBoundIntegrator
from gamd.langevin.dual_boost_integrators import LowerBoundIntegrator as DualLowerBoundIntegrator
from gamd.langevin.dual_boost_integrators import UpperBoundIntegrator as DualUpperBoundIntegrator
from gamd.integrator_factory import *


def load_pdb_positions_and_box_vectors(file_type,
                                       coords_filename, need_box):
    if file_type == "pdb":
        positions = openmm_app.PDBFile(coords_filename)
    elif file_type == "mmcif":
        positions = openmm_app.PDBxFile(coords_filename)
    else:
        message = "File Type " + str(file_type) + \
                  " not implemented for loading positions and vox vectors."
        raise NotImplementedError(message)
    
    pdb_parmed = parmed.load_file(coords_filename)
    if need_box:
        assert pdb_parmed.box_vectors is not None, "No box vectors "\
            "found in {}. ".format(coords_filename) \
            + "Box vectors for an anchor must be defined with a CRYST "\
            "line within the PDB file."

    return positions, pdb_parmed.box_vectors


class GamdSimulation:
    def __init__(self):
        self.system = None
        self.integrator = None
        self.simulation = None
        self.traj_reporter = None
        self.first_boost_group = None
        self.second_boost_group = None
        self.first_boost_type = None
        self.second_boost_type = None
        self.platform = "CUDA"
        self.device_index = 0


class GamdSimulationFactory:
    def __init__(self):
        return

    @staticmethod
    def configure_platform(gamd_simulation, topology,
                           platform_name, device_index):
        properties = {}
        user_platform_name = platform_name.lower()
        #
        # NOTE:  The Platform names are case sensitive.  From the OpenMM
        # docs "The platform name should be one of OpenCL, CUDA, CPU,
        # or Reference."
        #
        if user_platform_name == "cuda":
            platform = openmm.Platform.getPlatformByName('CUDA')
            properties['CudaPrecision'] = 'mixed'
            properties['DeviceIndex'] = device_index
            gamd_simulation.simulation = openmm_app.Simulation(
                topology.topology, gamd_simulation.system,
                gamd_simulation.integrator, platform, properties)
            gamd_simulation.device_index = device_index
            gamd_simulation.platform = 'CUDA'
        elif user_platform_name == "opencl":
            platform = openmm.Platform.getPlatformByName('OpenCL')
            properties['DeviceIndex'] = device_index
            gamd_simulation.simulation = openmm_app.Simulation(
                topology.topology, gamd_simulation.system,
                gamd_simulation.integrator, platform, properties)
            gamd_simulation.device_index = device_index
            gamd_simulation.platform = 'OpenCL'
        else:
            platform = openmm.Platform.getPlatformByName(platform_name)
            gamd_simulation.simulation = openmm_app.Simulation(
                topology.topology, gamd_simulation.system,
                gamd_simulation.integrator, platform)
            gamd_simulation.platform = platform_name

    def configure_barostat(config, gamd_simulation):
        if config.barostat is not None:
            barostat = openmm.MonteCarloBarostat(
                config.barostat.pressure,
                config.temperature,
                config.barostat.frequency)
            gamd_simulation.system.addForce(barostat)


    @staticmethod
    def configure_gamd_langevin_integrator(config, gamd_simulation):
        boost_type_str = config.integrator.boost_type
        gamdIntegratorFactory = GamdIntegratorFactory()
        result = gamdIntegratorFactory.get_integrator(config, system)

        [gamd_simulation.first_boost_group,
         gamd_simulation.second_boost_group,
         integrator, gamd_simulation.first_boost_type,
         gamd_simulation.second_boost_type] = result

        integrator.setRandomNumberSeed(config.integrator.random_seed)
        integrator.setFriction(config.integrator.friction_coefficient)
        gamd_simulation.integrator = integrator


    @staticmethod
    def create_amber_simulation(config, need_box,
                                nonbonded_method,
                                constraints):
        box_vectors = None
        gamd_simulation = GamdSimulation()
        prmtop = openmm_app.AmberPrmtopFile(
            config.input_files.amber.topology)
        topology = prmtop
        if config.input_files.amber.coordinates_filetype in ["inpcrd",
                                                             "rst7"]:
            positions = openmm_app.AmberInpcrdFile(
                config.input_files.amber.coordinates)
            box_vectors = positions.boxVectors
        elif config.input_files.amber.coordinates_filetype in ["pdb", "mmcif"]:
            coords_filename = config.input_files.amber.coordinates
            file_type = config.input_files.amber.coordinates_filetype
            positions, \
                box_vectors = load_pdb_positions_and_box_vectors(file_type,
                                                                 coords_filename,
                                                                 need_box)
        else:
            raise Exception("Invalid input type: %s. Allowed types are: " \
                            "'pdb' and 'rst7'/'inpcrd'.")
        gamd_simulation.system = prmtop.createSystem(
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=config.system.nonbonded_cutoff,
            constraints=constraints)

        return gamd_simulation, topology, positions, box_vectors

    @staticmethod
    def create_charmm_simulation(config,
                                 nonbonded_method,
                                 constraints):
        box_vecotrs = None
        gamd_simulation = GamdSimulation()
        psf = openmm_app.CharmmPsfFile(config.input_files.charmm.topology)

        if config.input_files.charmm.coordinates_filetype == "crd":
            positions = openmm_app.CharmmCrdFile(
                config.input_files.charmm.coordinates)
        elif config.input_files.charmm.coordinates_filetype == "pdb":
            positions = openmm_app.PDBFile(
                config.input_files.charmm.coordinates)
        elif config.input_files.charmm.coordinates_filetype == "mmcif":
            positions = openmm_app.PDBxFile(
                config.input_files.charmm.coordinates)
        else:
            raise Exception("Invalid input type: %s. Allowed types are: " \
                            "'crd' and 'pdb'.")

        # if a custom set of box vectors were defined in the configuration file,
        # then we need to set it in the psf object prior to system creation
        # from the psf object.
        if config.input_files.charmm.is_config_box_vector_defined:
            psf.setBox(*config.input_files.charmm.box_vectors)

        params = openmm_app.CharmmParameterSet(
            *config.input_files.charmm.parameters)

        topology = psf
        gamd_simulation.system = psf.createSystem(
            params=params,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=config.system.nonbonded_cutoff,
            switchDistance=config.system.switch_distance,
            ewaldErrorTolerance=config.system.ewald_error_tolerance,
            constraints=constraints)

        return gamd_simulation, topology, positions, box_vectors

    @staticmethod
    def create_gromacs_simulation(config, nonbonded_method, constraints):
        gamd_simulation = GamdSimulation()
        gro = openmm_app.GromacsGroFile(
            config.input_files.gromacs.coordinates)
        top = openmm_app.GromacsTopFile(
            config.input_files.gromacs.topology,
            periodicBoxVectors=gro.getPeriodicBoxVectors(),
            includeDir=config.input_files.gromacs.include_dir)
        box_vectors = gro.getPeriodicBoxVectors()
        topology = top
        positions = gro
        gamd_simulation.system = top.createSystem(
            nonbondedMethod=nonbonded_method,
            nonbondedCutoff=config.system.nonbonded_cutoff,
            constraints=constraints)

        return gamd_simulation, topology, positions, box_vectors

    @staticmethod
    def create_forcefield_simulation(config, need_box,
                                     nonbonded_method,
                                     constraints):
        gamd_simulation = GamdSimulation()
        coords_filename = config.input_files.forcefield.coordinates
        file_type = config.input_files.forcefield.coordinates_filetype

        positions, \
            box_vectors = load_pdb_positions_and_box_vectors(file_type,
                                                             coords_filename,
                                                             need_box)

        forcefield_filenames \
            = config.input_files.forcefield.forcefield_list_native \
              + config.input_files.forcefield.forcefield_list_external

        forcefield = openmm_app.ForceField(*forcefield_filenames)
        topology = positions
        gamd_simulation.system = forcefield.createSystem(
            topology.topology,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=config.system.nonbonded_cutoff,
            constraints=constraints)

        return gamd_simulation, topology, positions, box_vectors


    @staticmethod
    def get_nonbonded_method(config):
        result = None
        need_box = True
        method_str = config.system.nonbonded_method
        if method_str == "pme":
            result = openmm_app.PME
        elif method_str == "nocutoff":
            result = openmm_app.NoCutoff
            need_box = False
        elif method_str == "cutoffnonperiodic":
            result = openmm_app.CutoffNonPeriodic
        elif method_str == "cutoffperiodic":
            result = openmm_app.CutoffPeriodic
        elif method_str == "ewald":
            result = openmm_app.Ewald
        else:
            raise Exception("nonbonded method not found: %s",
                            method_str)

        return result, need_box

    @staticmethod
    def get_constraints(config):
        result = None
        constraints_str = config.system.constraints
        if constraints_str == "none" or constraints_str is None:
            result = None
        elif constraints_str == "hbonds":
            result = openmm_app.HBonds
        elif constraints_str == "allbonds":
            result = openmm_app.AllBonds
        elif constraints_str == "hangles":
            result = openmm_app.HAngles
        else:
            raise Exception("constraints not found: %s",
                            constraints_str)

        return result


    def _handle_inputs(self, config, need_box, nonbonded_method, constraints):
        gamd_simulation = None
        topology = None
        positions = None
        box_vectors = None
        if config.input_files.amber is not None:
            gamd_simulation, \
                topology, \
                positions, \
                box_vectors = self.create_amber_simulation(config,
                                                           need_box,
                                                           nonbonded_method,
                                                           constraints)
        elif config.input_files.charmm is not None:
            gamd_simulation, \
                topology, \
                positions, \
                box_vectors = self.create_charmm_simulation(config,
                                                            nonbonded_method,
                                                            constraints)

        elif config.input_files.gromacs is not None:
            gamd_simulation, \
                topology, \
                positions, \
                box_vectors = self.create_gromacs_simulation(config,
                                                             nonbonded_method,
                                                             constraints)

        elif config.input_files.forcefield is not None:
            gamd_simulation, \
                topology, \
                positions, \
                box_vectors = self.create_forcefield_simulation(config,
                                                           need_box,
                                                           nonbonded_method,
                                                           constraints)
        else:
            raise Exception("No valid input files found. OpenMM simulation "
                            "not made.")

        return gamd_simulation, topology, positions, box_vectors

    @staticmethod
    def configure_trajectory_reporting(config, gamd_simulation):
        if config.outputs.reporting.coordinates_file_type == "dcd":
            gamd_simulation.traj_reporter = openmm_app.DCDReporter

        elif config.outputs.reporting.coordinates_file_type == "pdb":
            gamd_simulation.traj_reporter = openmm_app.PDBReporter

        else:
            raise Exception("Reporter type not found:",
                            config.outputs.reporting.coordinates_file_type)

    def createGamdSimulation(self, config, platform_name, device_index):
        nonbonded_method, need_box = self.get_nonbonded_method(config)
        constraints = self.get_constraints(config)

        gamd_simulation, \
            topology, \
            positions, \
            box_vectors = self._handle_inputs(config,
                                              need_box,
                                              nonbonded_method,
                                              constraints)

        if config.integrator.algorithm == "langevin":
            self.configure_gamd_langevin_integrator(config,
                                                    gamd_simulation)
        else:
            raise Exception("Algorithm not implemented:",
                            config.integrator.algorithm)

        self.configure_barostat(config, gamd_simulation)
        self.configure_platform(gamd_simulation, topology,
                                platform_name, device_index)

        simulation = gamd_simulation.simulation
        simulation.context.setPositions(positions.positions)

        #
        # If this isn't a charmm configuration, but box vectors were defined,
        # then we can setup the box vectors after the simulation
        # object has been created.  (charmm psf requires it prior to the
        # simulation creation.)
        #
        if box_vectors is not None and config.input_files.charmm is None:
            gamd_simulation.simulation.context.setPeriodicBoxVectors(
                *box_vectors)

        if config.run_minimization:
            gamd_simulation.simulation.minimizeEnergy()

        simulation.context.setVelocitiesToTemperature(config.temperature)
        self.configure_trajectory_reporting(config, gamd_simulation)
    
        return gamd_simulation


if __name__ == "__main__":
    pass
