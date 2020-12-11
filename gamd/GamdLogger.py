from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from .stage_integrator import BoostType
from .stage_integrator import GamdStageIntegrator


class GamdLogger:

    def __init__(self, filename, mode, integrator, simulation):
        self.filename = filename
        self.gamdLog = open(filename, mode)
        self.integrator = integrator
        self.simulation = simulation

        self.starting_total_potential_energy = 0
        self.starting_dihedral_potential_energy = 0

        self.total_force_scaling_factor_name = integrator.get_variable_name_by_type(BoostType.TOTAL,
                                                                                    "ForceScalingFactor")
        self.dihedral_force_scaling_factor_name = integrator.get_variable_name_by_type(BoostType.DIHEDRAL,
                                                                                       "ForceScalingFactor")
        self.total_boost_potential_name = integrator.get_variable_name_by_type(BoostType.TOTAL, "BoostPotential")
        self.dihedral_boost_potential_name = integrator.get_variable_name_by_type(BoostType.DIHEDRAL, "BoostPotential")

    def __del__(self):
        self.gamdLog.close()

    def write_header(self):
        self.gamdLog.write("# Gaussian accelerated Molecular Dynamics log file\n")
        self.gamdLog.write("# All energy terms are stored in unit of kcal/mol\n")
        self.gamdLog.write(
            "# ntwx,total_nstep,Unboosted-Potential-Energy,Unboosted-Dihedral-Energy,Total-Force-Weight,Dihedral-Force-Weight,Boost-Energy-Potential,Boost-Energy-Dihedral\n")

    def mark_energies(self, group=None):
        state = self.simulation.context.getState(getEnergy=True)
        self.starting_total_potential_energy = state.getPotentialEnergy()
        if group is not None:
            dihedral_state = self.simulation.context.getState(getEnergy=True, groups={group})
            self.starting_dihedral_potential_energy = dihedral_state.getPotentialEnergy()
        else:
            self.starting_dihedral_potential_energy = 0 * kilojoules_per_mole

    def write_to_gamd_log(self, step):
        force_scaling_factors = self.integrator.get_force_scaling_factors()
        boost_potentials = self.integrator.get_boost_potentials()

        total_potential_energy = str(self.starting_total_potential_energy / (kilojoules_per_mole * 4.184))
        dihedral_energy = str(self.starting_dihedral_potential_energy / (kilojoules_per_mole * 4.184))
        total_force_scaling_factor = str(force_scaling_factors[self.total_force_scaling_factor_name])
        dihedral_force_scaling_factor = str(force_scaling_factors[self.dihedral_force_scaling_factor_name])
        total_boost_potential = str(boost_potentials[self.total_boost_potential_name] / 4.184)
        dihedral_boost_potential = str(boost_potentials[self.dihedral_boost_potential_name] / 4.184)

        self.gamdLog.write("\t" + str(1) + "\t" + str(step * 1) + "\t" +
                           total_potential_energy + "\t" +
                           dihedral_energy + "\t" +
                           total_force_scaling_factor + "\t" +
                           dihedral_force_scaling_factor + "\t" +
                           total_boost_potential + "\t" +
                           dihedral_boost_potential + "\n")
