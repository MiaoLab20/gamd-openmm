'''
Created on Oct 28, 2020

Creates all OpenMM objects from Config() object that can be used in a 
GaMD simulation.

@author: lvotapka
'''
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


def load_pdb_positions_and_box_vectors(pdb_coords_filename):
    positions = openmm_app.PDBFile(pdb_coords_filename)
    pdb_parmed = parmed.load_file(pdb_coords_filename)
    assert pdb_parmed.box_vectors is not None, "No box vectors "\
        "found in {}. ".format(pdb_coords_filename) \
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
        self.cuda_device_index = "0"


class GamdSimulationFactory:
    def __init__(self):
        return
        
    def createGamdSimulation(self, config, platform_name, device_index):
        if config.system.nonbonded_method == "pme":
            nonbondedMethod = openmm_app.PME
            
        elif config.system.nonbonded_method == "nocutoff":
            nonbondedMethod = openmm_app.NoCutoff
            
        elif config.system.nonbonded_method == "cutoffnonperiodic":
            nonbondedMethod = openmm_app.CutoffNonPeriodic
            
        elif config.system.nonbonded_method == "cutoffperiodic":
            nonbondedMethod = openmm_app.CutoffPeriodic
            
        elif config.system.nonbonded_method == "ewald":
            nonbondedMethod = openmm_app.Ewald
        
        else:
            raise Exception("nonbonded method not found: %s", 
                            config.system.nonbonded_method)
        
        if config.system.constraints == "none" \
                or config.system.constraints is None:
            constraints = None
        
        elif config.system.constraints == "hbonds":
            constraints = openmm_app.HBonds
            
        elif config.system.constraints == "allbonds":
            constraints = openmm_app.AllBonds
            
        elif config.system.constraints == "hangles":
            constraints = openmm_app.HAngles
            
        else:
            raise Exception("constraints not found: %s", 
                            config.system.constraints)
            
        box_vectors = None
        gamdSimulation = GamdSimulation()
        if config.input_files.amber is not None:
            prmtop = openmm_app.AmberPrmtopFile(
                config.input_files.amber.topology)
            topology = prmtop
            if config.input_files.amber.coordinates_filetype in ["inpcrd", 
                                                                 "rst7"]:
                positions = openmm_app.AmberInpcrdFile(
                    config.input_files.amber.coordinates)
                box_vectors = positions.boxVectors
            elif config.input_files.amber.coordinates_filetype == "pdb":
                pdb_coords_filename = config.input_files.amber.coordinates
                positions, box_vectors = load_pdb_positions_and_box_vectors(
                    pdb_coords_filename)
            else:
                raise Exception("Invalid input type: %s. Allowed types are: "\
                                "'pdb' and 'rst7'/'inpcrd'.")
            gamdSimulation.system = prmtop.createSystem(
                nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=config.system.nonbonded_cutoff, 
                constraints=constraints)
            
        elif config.input_files.charmm is not None:
            psf = openmm_app.CharmmPsfFile(
                config.input_files.charmm.topology)
            pdb_coords_filename = config.input_files.charmm.coordinates
            positions, box_vectors = load_pdb_positions_and_box_vectors(
                pdb_coords_filename)
            params = openmm_app.CharmmParameterSet(
                config.input_files.charmm.parameters)
            topology = psf
            gamdSimulation.system = psf.createSystem(
                params=params,
                nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=config.system.nonbonded_cutoff, 
                constraints=constraints)
            
        elif config.input_files.gromacs is not None:
            gro = openmm_app.GromacsGroFile(
                config.input_files.gromacs.coordinates)
            top = openmm_app.GromacsTopFile(
                config.input_files.gromacs.topology,
                periodicBoxVectors=gro.getPeriodicBoxVectors(),
                includeDir=config.input_files.gromacs.include_dir)
            topology = top
            positions = gro
            gamdSimulation.system = top.createSystem(
                nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=config.system.nonbonded_cutoff, 
                constraints=constraints)
            
        elif config.input_files.forcefield is not None:
            pdb_coords_filename = config.input_files.forcefield.coordinates
            positions, box_vectors = load_pdb_positions_and_box_vectors(
                pdb_coords_filename)
            forcefield_filenames \
                = config.input_files.forcefield.forcefield_list_native \
                + config.input_files.forcefield.forcefield_list_external
            forcefield = openmm_app.ForceField(*forcefield_filenames)
            topology = positions
            gamdSimulation.system = forcefield.createSystem(
                topology.topology,
                nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=config.system.nonbonded_cutoff, 
                constraints=constraints)
            
        else:
            raise Exception("No valid input files found. OpenMM simulation "\
                            "not made.")
        if config.integrator.algorithm == "langevin":
            boost_type_str = config.integrator.boost_type
            gamdIntegratorFactory = GamdIntegratorFactory()
            result = gamdIntegratorFactory.get_integrator(
                boost_type_str, gamdSimulation.system, config.temperature, 
                config.integrator.dt, 
                config.integrator.number_of_steps.conventional_md_prep, 
                config.integrator.number_of_steps.conventional_md, 
                config.integrator.number_of_steps.gamd_equilibration_prep, 
                config.integrator.number_of_steps.gamd_equilibration, 
                config.integrator.number_of_steps.total_simulation_length, 
                config.integrator.number_of_steps.averaging_window_interval,
                sigma0p=config.integrator.sigma0.primary, 
                sigma0d=config.integrator.sigma0.secondary)
            [gamdSimulation.first_boost_group, 
             gamdSimulation.second_boost_group, 
             integrator, gamdSimulation.first_boost_type, 
             gamdSimulation.second_boost_type] = result
            integrator.setRandomNumberSeed(config.integrator.random_seed)
            integrator.setFriction(config.integrator.friction_coefficient)
            gamdSimulation.integrator = integrator
           
        else:
            raise Exception("Algorithm not implemented:", 
                            config.integrator.algorithm)
        
        if config.barostat is not None:
            barostat = openmm.MonteCarloBarostat(
                config.barostat.pressure, 
                config.temperature,
                config.barostat.frequency)
            gamdSimulation.system.addForce(barostat)
        
        properties = {}
        if platform_name == "cuda":
            platform = openmm.Platform.getPlatformByName(platform_name)
            properties['CudaPrecision'] = "mixed"
            properties['DeviceIndex'] = device_index
            gamdSimulation.simulation = openmm_app.Simulation(
                topology.topology, gamdSimulation.system, 
                gamdSimulation.integrator, platform, properties)
        elif platform_name == "opencl":
            platform = openmm.Platform.getPlatformByName(platform_name)
            properties['DeviceIndex'] = device_index
            gamdSimulation.simulation = openmm_app.Simulation(
                topology.topology, gamdSimulation.system, 
                gamdSimulation.integrator, platform, properties)
        else:
            gamdSimulation.simulation = openmm_app.Simulation(
                topology.topology, gamdSimulation.system, 
                gamdSimulation.integrator)
        
        gamdSimulation.simulation.context.setPositions(positions.positions)
        gamdSimulation.simulation.context.setPeriodicBoxVectors(
            *box_vectors)
        if config.run_minimization:
            gamdSimulation.simulation.minimizeEnergy()
        
        gamdSimulation.simulation.context.setVelocitiesToTemperature(
            config.temperature)
        
        if config.outputs.reporting.coordinates_file_type == "dcd":
            gamdSimulation.traj_reporter = openmm_app.DCDReporter
            
        elif config.outputs.reporting.coordinates_file_type == "pdb":
            gamdSimulation.traj_reporter = openmm_app.PDBReporter
            
        else:
            raise Exception("Reporter type not found:", 
                            config.outputs.reporting.coordinates.file_type)
    
        return gamdSimulation

    
if __name__ == "__main__":
    pass