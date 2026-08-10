"""
Created on Oct 28, 2020

Creates all OpenMM objects from Config() object that can be used in a
GaMD simulation.

Original author:
    lvotapka

Modified by:
    Alaa REGUEI (2026)

Modifications:
    - Adapted for membrane protein GPCR simulations
    - Added Hydrogen Mass Repartitioning (HMR) support to enable 4 fs integration timestep
    - Added semi-isotropic barostat conditions for membrane systems (monte carlo membrane barostat)
    - Ability to initialize simulations from an OpenMM state XML file
      (e.g., previously equilibrated systems or systems prepared with OpenMM)
    - reporting of additional simulation data such as potential energy, total energy, temperature, and simulation speed (ns/day) at user-specified intervals
"""

import os
import sys

import parmed
import openmm as openmm
import openmm.app as openmm_app
import openmm.unit as unit
from openmm import XmlSerializer

from gamd import parser
from gamd.langevin.total_boost_integrators import (
    LowerBoundIntegrator as TotalLowerBoundIntegrator,
)
from gamd.langevin.total_boost_integrators import (
    UpperBoundIntegrator as TotalUpperBoundIntegrator,
)
from gamd.langevin.dihedral_boost_integrators import (
    LowerBoundIntegrator as DihedralLowerBoundIntegrator,
)
from gamd.langevin.dihedral_boost_integrators import (
    UpperBoundIntegrator as DihedralUpperBoundIntegrator,
)
from gamd.langevin.dual_boost_integrators import (
    LowerBoundIntegrator as DualLowerBoundIntegrator,
)
from gamd.langevin.dual_boost_integrators import (
    UpperBoundIntegrator as DualUpperBoundIntegrator,
)
from gamd.integrator_factory import *


def load_pdb_positions_and_box_vectors(pdb_coords_filename, need_box):
    pdb = openmm_app.PDBFile(pdb_coords_filename)
    pdb_parmed = parmed.load_file(pdb_coords_filename)

    if need_box:
        assert pdb_parmed.box_vectors is not None, (
            f"No box vectors found in {pdb_coords_filename}. "
            "Box vectors for an anchor must be defined with a CRYST line within the PDB file."
        )

    return pdb.positions, pdb_parmed.box_vectors


def load_xml_positions_and_box_vectors(xml_filename, need_box=True):
    with open(xml_filename, "r") as f:
        state = XmlSerializer.deserialize(f.read())

    positions = state.getPositions()

    box_vectors = None
    if need_box:
        box_vectors = state.getPeriodicBoxVectors()
        assert (
            box_vectors is not None
        ), f"No periodic box vectors found in {xml_filename}."
    return positions, box_vectors


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
        self.topology = None
        self.positions = None
        self.box_vectors = None
        self.platform = "CUDA"
        self.device_index = 0


class GamdSimulationFactory:
    def __init__(self):
        return

    def createGamdSimulation(self, config, platform_name, device_index):

        # -------------------------
        # NONBONDED METHOD
        # -------------------------
        need_box = True
        if config.system.nonbonded_method == "pme":
            nonbondedMethod = openmm_app.PME
        elif config.system.nonbonded_method == "nocutoff":
            nonbondedMethod = openmm_app.NoCutoff
            need_box = False
        elif config.system.nonbonded_method == "cutoffperiodic":
            nonbondedMethod = openmm_app.CutoffPeriodic
        elif config.system.nonbonded_method == "ewald":
            nonbondedMethod = openmm_app.Ewald
        else:
            raise Exception("nonbonded method not found")

        # -------------------------
        # CONSTRAINTS
        # -------------------------
        if config.system.constraints in [None, "none"]:
            constraints = None
        elif config.system.constraints == "hbonds":
            constraints = openmm_app.HBonds
        elif config.system.constraints == "allbonds":
            constraints = openmm_app.AllBonds
        else:
            raise Exception("constraints not found")

        # -------------------------
        # LOAD TOPOLOGY AND POSITIONS
        # -------------------------
        prmtop = openmm_app.AmberPrmtopFile(config.input_files.amber.topology)
        if config.input_files.amber.coordinates_filetype in ["inpcrd", "rst7"]:
            positions = openmm_app.AmberInpcrdFile(config.input_files.amber.coordinates)
            box_vectors = positions.boxVectors
        elif config.input_files.amber.coordinates_filetype == "pdb":
            pdb_coords_filename = config.input_files.amber.coordinates
            positions, box_vectors = load_pdb_positions_and_box_vectors(
                pdb_coords_filename, need_box
            )
        elif config.input_files.amber.coordinates_filetype == "xml":
            xml_coords_filename = config.input_files.amber.coordinates
            positions, box_vectors = load_xml_positions_and_box_vectors(
                xml_coords_filename, need_box
            )

        # -------------------------
        # SYSTEM CREATION with HMR (if specified)
        # -------------------------
        system = prmtop.createSystem(
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=config.system.nonbonded_cutoff,
            constraints=constraints,
            rigidWater=True,
            hydrogenMass=1.5 * unit.amu,
            ewaldErrorTolerance=0.0005,
        )

        if box_vectors is not None:
            print("Box vectors:", box_vectors)
            system.setDefaultPeriodicBoxVectors(*box_vectors)

        gamdSimulation = GamdSimulation()
        gamdSimulation.system = system
        gamdSimulation.topology = prmtop.topology
        gamdSimulation.positions = positions
        gamdSimulation.box_vectors = box_vectors

        # -------------------------
        # INTEGRATOR
        # -------------------------
        if config.integrator.algorithm != "langevin":
            raise Exception("Algorithm not implemented")

        gamdIntegratorFactory = GamdIntegratorFactory()
        result = gamdIntegratorFactory.get_integrator(
            config.integrator.boost_type,
            system,
            config.temperature,
            config.integrator.dt,
            config.integrator.number_of_steps.conventional_md_prep,
            config.integrator.number_of_steps.conventional_md,
            config.integrator.number_of_steps.gamd_equilibration_prep,
            config.integrator.number_of_steps.gamd_equilibration,
            config.integrator.number_of_steps.total_simulation_length,
            config.integrator.number_of_steps.averaging_window_interval,
            sigma0p=config.integrator.sigma0.primary,
            sigma0d=config.integrator.sigma0.secondary,
        )

        (
            gamdSimulation.first_boost_group,
            gamdSimulation.second_boost_group,
            integrator,
            gamdSimulation.first_boost_type,
            gamdSimulation.second_boost_type,
        ) = result

        integrator.setRandomNumberSeed(config.integrator.random_seed)
        integrator.setFriction(config.integrator.friction_coefficient)

        gamdSimulation.integrator = integrator

        # -------------------------
        # BAROSTAT
        # -------------------------
        if config.barostat is not None:
            barostat = openmm.MonteCarloMembraneBarostat(
                1 * unit.bar,
                0 * unit.bar * unit.nanometer,
                config.temperature,
                openmm.MonteCarloMembraneBarostat.XYIsotropic,
                openmm.MonteCarloMembraneBarostat.ZFree,
                100,
            )
            system.addForce(barostat)

        # -------------------------
        # PLATFORM
        # -------------------------
        platform = openmm.Platform.getPlatformByName(platform_name)

        gamdSimulation.simulation = openmm_app.Simulation(
            gamdSimulation.topology, system, integrator, platform
        )
        gamdSimulation.simulation.context.setPositions(gamdSimulation.positions)

        gamdSimulation.simulation.context.setPeriodicBoxVectors(*box_vectors)
        print(
            "Initial potential energy: ",
            gamdSimulation.simulation.context.getState(
                getEnergy=True
            ).getPotentialEnergy(),
        )

        if config.run_minimization:
            print("Running energy minimization...")
            gamdSimulation.simulation.minimizeEnergy(maxIterations=10000)
        print(
            "Potential energy after minimization: ",
            gamdSimulation.simulation.context.getState(
                getEnergy=True
            ).getPotentialEnergy(),
        )

        gamdSimulation.simulation.context.setVelocitiesToTemperature(config.temperature)

        # -------------------------
        # REPORTERS
        # -------------------------
        if config.outputs.reporting.coordinates_file_type == "dcd":
            gamdSimulation.traj_reporter = openmm_app.DCDReporter
        elif config.outputs.reporting.coordinates_file_type == "pdb":
            gamdSimulation.traj_reporter = openmm_app.PDBReporter
        else:
            raise Exception("Reporter type not found")
            
        report_interval = getattr(config.outputs.reporting, "energy_interval", 1000)
        
        reporter = openmm_app.StateDataReporter(
            sys.stdout,             # File or stream to write to
            report_interval,         # Interval in steps
            step=True,               # Prints current step
            potentialEnergy=True,    # Prints potential energy
            totalEnergy=True,
            temperature=True,        # Prints temperature
            speed=True,              # Calculates and prints simulation speed (ns/day)
            separator='\t',
        )
        
        # Attach the reporter directly to the OpenMM Simulation object
        gamdSimulation.simulation.reporters.append(reporter)        
        return gamdSimulation

if __name__ == "__main__":
    pass
