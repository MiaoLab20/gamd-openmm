'''
Created on Oct 28, 2020

Creates all OpenMM objects from Config() object that can be used in a 
GaMD simulation.

@author: lvotapka
'''
import os

from simtk import openmm
import simtk.openmm.app as openmm_app
from simtk import unit

from gamd import parser
# change to generic integrator someday
from gamd.langevin.total_boost_integrators import LowerBoundIntegrator as TotalLowerBoundIntegrator
from gamd.langevin.total_boost_integrators import UpperBoundIntegrator as TotalUpperBoundIntegrator
from gamd.langevin.dihedral_boost_integrators import LowerBoundIntegrator as DihedralLowerBoundIntegrator
from gamd.langevin.dihedral_boost_integrators import UpperBoundIntegrator as DihedralUpperBoundIntegrator

class GamdSimulation:
    def __init__(self):
        self.system = None
        self.integrator = None
        self.simulation = None
        self.traj_reporter = None
        
    
class GamdSimulationFactory:
    def __init__(self):
        return
        
    def createGamdSimulation(self, config):
        gamdSimulation = GamdSimulation()
        if config.system_files_config.type == "amber":
            prmtop = openmm_app.AmberPrmtopFile(
                config.system_files_config.prmtop_filename)
            inpcrd = openmm_app.AmberInpcrdFile(
                config.system_files_config.inpcrd_filename)
            if config.system_files_config.load_box_vecs_from_coords_file:
                config.box_vectors = inpcrd.boxVectors
            topology = prmtop
            positions = inpcrd
            
        elif config.system_files_config.type == "charmm":
            psf = openmm_app.CharmmPsfFile(
                config.system_files_config.psf_filename)
            pdb = openmm_app.PDBFile(
                config.system_files_config.pdb_filename)
            params = openmm_app.CharmmParameterSet(
                config.system_files_config.params_filenames)
            topology = psf
            positions = pdb
            
        elif config.system_files_config.type == "gromacs":
            gro = openmm_app.GromacsGroFile(
                config.system_files_config.gro_filename)
            top = openmm_app.GromacsTopFile(
                config.system_files_config.top_filename,
                periodicBoxVectors=gro.getPeriodicBoxVectors(),
                includeDir=config.system_files_config.include_dir)
            topology = top
            positions = gro
            
        elif config.system_files_config.type == "forcefield":
            pdb = openmm_app.PDBFile(
                config.system_files_config.pdb_filename)
            forcefield = openmm_app.ForceField(
                config.system_files_config.forcefield_filenames)
            topology = pdb
            positions = pdb
        
        else:
            print("Type:", config.system_files_config.type, "not found.",
                  "OpenMM files not made.")
        
        if config.nonbonded_method == "pme":
            nonbondedMethod = openmm_app.PME
            
        elif config.nonbonded_method == "nocutoff":
            nonbondedMethod = openmm_app.NoCutoff
            
        elif config.nonbonded_method == "cutoffnonperiodic":
            nonbondedMethod = openmm_app.CutoffNonPeriodic
            
        elif config.nonbonded_method == "cutoffperiodic":
            nonbondedMethod = openmm_app.CutoffPeriodic
            
        elif config.nonbonded_method == "ewald":
            nonbondedMethod = openmm_app.Ewald
        
        else:
            raise Exception("nonbonded method not found: %s", 
                            config.nonbonded_method)
        
        if config.constraints == "none" or config.constraints is None:
            constraints = None
        
        elif config.constraints == "hbonds":
            constraints = openmm_app.HBonds
            
        elif config.constraints == "allbonds":
            constraints = openmm_app.AllBonds
            
        elif config.constraints == "hangles":
            constraints = openmm_app.HAngles
            
        else:
            raise Exception("constraints not found: %s", 
                            config.constraints)
        
        if config.system_files_config.type == "amber":
            gamdSimulation.system = prmtop.createSystem(
                nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=config.nonbonded_cutoff, 
                constraints=constraints)
        
        elif config.system_files_config.type == "charmm":
            gamdSimulation.system = psf.createSystem(
                params=params,
                nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=config.nonbonded_cutoff, 
                constraints=constraints)
            
        elif config.system_files_config.type == "gromacs":
            gamdSimulation.system = top.createSystem(
                nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=config.nonbonded_cutoff, 
                constraints=constraints)
            
        elif config.system_files_config.type == "forcefield":
            gamdSimulation.system = forcefield.createSystem(
                pdb.topology,
                nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=config.nonbonded_cutoff, 
                constraints=constraints)
        
        else:
            print("Type:", config.system_files_config.type, "not found.",
                  "OpenMM files not made.")
        
        if config.integrator_type == "langevin":
            if config.total_boost and not config.dihedral_boost:
                if config.gamd_bound == 'lower':
                    gamdSimulation.integrator = TotalLowerBoundIntegrator(
                        dt=config.dt, 
                        ntcmdprep=config.num_steps_conventional_md_prep, 
                        ntcmd=config.num_steps_conventional_md, 
                        ntebprep=config.num_steps_gamd_equilibration_prep, 
                        nteb=config.num_steps_gamd_equilibration, 
                        nstlim=config.total_simulation_length, 
                        ntave=config.num_steps_per_averaging,
                        sigma0=config.total_boost_sigma0,
                        collision_rate=config.friction_coefficient,
                        temperature=config.target_temperature,
                        restart_filename=None)
                elif config.gamd_bound == 'upper':
                    gamdSimulation.integrator = TotalUpperBoundIntegrator(
                        dt=config.dt, 
                        ntcmdprep=config.num_steps_conventional_md_prep, 
                        ntcmd=config.num_steps_conventional_md, 
                        ntebprep=config.num_steps_gamd_equilibration_prep, 
                        nteb=config.num_steps_gamd_equilibration, 
                        nstlim=config.total_simulation_length, 
                        ntave=config.num_steps_per_averaging,
                        sigma0=config.total_boost_sigma0,
                        collision_rate=config.friction_coefficient,
                        temperature=config.target_temperature,
                        restart_filename=None)
                else:
                    raise Exception(
                        "This type of langevin integrator not implemented:",
                        config.gamd_bound)
            elif not config.total_boost and config.dihedral_boost:
                if config.gamd_bound == 'lower':
                    gamdSimulation.integrator = DihedralLowerBoundIntegrator(
                        group=config.dihedral_group,
                        dt=config.dt, 
                        ntcmdprep=config.num_steps_conventional_md_prep, 
                        ntcmd=config.num_steps_conventional_md, 
                        ntebprep=config.num_steps_gamd_equilibration_prep, 
                        nteb=config.num_steps_gamd_equilibration, 
                        nstlim=config.total_simulation_length, 
                        ntave=config.num_steps_per_averaging,
                        sigma0=config.dihedral_boost_sigma0,
                        collision_rate=config.friction_coefficient,
                        temperature=config.target_temperature,
                        restart_filename=None)
                        
                elif config.gamd_bound == 'upper':
                    gamdSimulation.integrator = DihedralUpperBoundIntegrator(
                        group=config.dihedral_group,
                        dt=config.dt, 
                        ntcmdprep=config.num_steps_conventional_md_prep, 
                        ntcmd=config.num_steps_conventional_md, 
                        ntebprep=config.num_steps_gamd_equilibration_prep, 
                        nteb=config.num_steps_gamd_equilibration, 
                        nstlim=config.total_simulation_length, 
                        ntave=config.num_steps_per_averaging,
                        sigma0=config.dihedral_boost_sigma0,
                        collision_rate=config.friction_coefficient,
                        temperature=config.target_temperature,
                        restart_filename=None)
                else:
                    raise Exception(
                        "This type of langevin integrator not implemented:",
                        config.gamd_bound)
            else:
                raise Exception(
                        "This combination of total and dihedral boosts "\
                        "not allowed or not yet implemented: Total_boost:",
                        config.total_boost, "Dihedral_boost:", 
                        config.dihedral_boost)
                
        else:
            raise Exception("Integrator type not implemented:", 
                            config.integrator_type)
        
        if config.use_barostat:
            barostat = openmm.MonteCarloBarostat(
                config.barostat_target_pressure, 
                config.barostat_target_temperature,
                config.barostat_frequency)
            gamdSimulation.system.addForce(barostat)
        
        # TODO: only dihedrals being boosted at this time.
        group = config.dihedral_group
        for force in gamdSimulation.system.getForces():
        #     print(force.__class__.__name__)
            if force.__class__.__name__ == 'PeriodicTorsionForce':
                force.setForceGroup(group)
                break        
        
        gamdSimulation.simulation = openmm_app.Simulation(
            topology.topology, gamdSimulation.system, 
            gamdSimulation.integrator)
        
        gamdSimulation.simulation.context.setPositions(positions.positions)
        gamdSimulation.simulation.context.setPeriodicBoxVectors(
            *config.box_vectors)
        if config.run_minimization:
            gamdSimulation.simulation.minimizeEnergy()
        
        gamdSimulation.simulation.context.setVelocitiesToTemperature(
            config.initial_temperature)
        
        if config.coordinates_reporter_file_type == None:
            gamdSimulation.traj_reporter = None
        
        elif config.coordinates_reporter_file_type == "dcd":
            gamdSimulation.traj_reporter = openmm_app.DCDReporter
            
        elif config.coordinates_reporter_file_type == "pdb":
            gamdSimulation.traj_reporter = openmm_app.PDBReporter
            
        else:
            raise Exception("Reporter type not found:", 
                            config.coordinates_reporter_file_type)
    
        return gamdSimulation
    
if __name__ == "__main__":
    parserFactory = parser.ParserFactory()
    config = parserFactory.parse_file("/home/lvotapka/gamd/sample_input.xml", "xml")
    gamdSimulationFactory = GamdSimulationFactory()
    gamdSimulation = gamdSimulationFactory.createGamdSimulation(config)