import openmm.unit as unit
from .stage_integrator import BoostType
from .stage_integrator import GamdStageIntegrator


class TrackedValue:
    def __init__(self, boost_type, group, tracked_integrator, tracked_simulation):
        self.__boost_type = boost_type
        self.__group = group
        self.__integrator = tracked_integrator
        self.__simulation = tracked_simulation
        self.__starting_potential_energy = 0 * unit.kilojoules_per_mole
        self.__boost_potential_name = tracked_integrator.get_variable_name_by_type(
            self.__boost_type, "BoostPotential")
        self.__force_scaling_factor_name = tracked_integrator.get_variable_name_by_type(
            self.__boost_type, "ForceScalingFactor")
        self.__effective_harmonic_constant_name = tracked_integrator.get_variable_name_by_type(
            self.__boost_type, "k0")

    def mark_energy(self):
        if self.__boost_type == BoostType.TOTAL:
            state = self.__simulation.context.getState(getEnergy=True)
            self.__starting_potential_energy = state.getPotentialEnergy()
        else:
            state = self.__simulation.context.getState(getEnergy=True, groups={self.__group})
            self.__starting_potential_energy = state.getPotentialEnergy()

    def get_reporting_force_scaling_factor(self):
        scaling_factors = self.__integrator.get_force_scaling_factors()
        return str(scaling_factors[self.__force_scaling_factor_name])

    def get_reporting_boost_potential(self):
        boost_potentials = self.__integrator.get_boost_potentials()
        return str(boost_potentials[self.__boost_potential_name] / 4.184)

    def get_reporting_starting_energy(self):
        return str(self.__starting_potential_energy /
                   (unit.kilojoules_per_mole * 4.184))

    def get_boost_type(self):
        return self.__boost_type

    def get_reporting_effective_harmonic_constant(self):
        effective_harmonic_constants = self.__integrator.get_effective_harmonic_constants()
        return str(effective_harmonic_constants[self.__effective_harmonic_constant_name])


class GamdLogger:

    def __init__(self, filename, mode, integrator, simulation,
                 first_boost_type, first_boost_group,
                 second_boost_type, second_boost_group):
        """
        Parameters
        ----------
        :param filename:           The gamd.log file path and file name.
        :param mode:               The write mode to output the file.
        :param integrator:         The integrator from which to pull information.
        :param simulation:         The simulation from which to pull information
        :param first_boost_type:   The simple boost type to record (no dual types)
        :param first_boost_group:  The group associated with the 1st boost type.  Empty double quoted string for total.
        :param second_boost_type:  The simple boost type to record (no dual types)
        :param second_boost_group: The group associated with the 2nd boost type.  Empty double quoted string for total.

        """

        self.filename = filename
        self.gamdLog = open(filename, mode)
        self.integrator = integrator
        self.simulation = simulation
        self.tracked_values = []

        if first_boost_type == BoostType.DUAL_TOTAL_DIHEDRAL or second_boost_type == BoostType.DUAL_TOTAL_DIHEDRAL:
            raise ValueError("The GamdLogger expects single value boost types as arguments, not compound boost types."
                             "  Compound boost types should be broken up.\n")

        if second_boost_type == BoostType.TOTAL:
            self.tracked_values.append(TrackedValue(second_boost_type, second_boost_group, integrator, simulation))
            self.tracked_values.append(TrackedValue(first_boost_type, first_boost_group, integrator, simulation))
        else:
            self.tracked_values.append(TrackedValue(first_boost_type, first_boost_group, integrator, simulation))
            self.tracked_values.append(TrackedValue(second_boost_type, second_boost_group, integrator, simulation))

    def __del__(self):
        self.gamdLog.close()

    def write_header(self):
        self.gamdLog.write("# Gaussian accelerated Molecular Dynamics log file\n")
        self.gamdLog.write("# All energy terms are stored in unit of kcal/mol\n")
        header_str = "# ntwx,total_nstep,Unboosted-{0}-Energy,Unboosted-{1}-Energy,{0}-Force-Weight,{1}-Force-Weight,{0}-Boost-Energy-Potential,{1}-Boost-Energy,{0}-Effctive-Harmonic-Constant,{1}-Effctive-Harmonic-Constant\n"
        header = header_str.format(self.tracked_values[0].get_boost_type().value,
                                   self.tracked_values[1].get_boost_type().value)
        self.gamdLog.write(header)

    def mark_energies(self):
        for tracked_value in self.tracked_values:
            tracked_value.mark_energy()

    def write_to_gamd_log(self, step):
        first_energy = self.tracked_values[0].get_reporting_starting_energy()
        second_energy = self.tracked_values[1].get_reporting_starting_energy()

        first_force_scaling_factor = self.tracked_values[0].get_reporting_force_scaling_factor()
        second_force_scaling_factor = self.tracked_values[1].get_reporting_force_scaling_factor()

        first_boost_potential = self.tracked_values[0].get_reporting_boost_potential()
        second_boost_potential = self.tracked_values[1].get_reporting_boost_potential()

        first_effective_harmonic_constant = self.tracked_values[0].get_reporting_effective_harmonic_constant()
        second_effective_harmonic_constant = self.tracked_values[1].get_reporting_effective_harmonic_constant()

        self.gamdLog.write("\t" + str(1) + "\t" + str(step * 1) + "\t" +
                           first_energy + "\t" +
                           second_energy + "\t" +
                           first_force_scaling_factor + "\t" +
                           second_force_scaling_factor + "\t" +
                           first_boost_potential + "\t" +
                           second_boost_potential + "\t" +
                           first_effective_harmonic_constant + "\t" +
                           second_effective_harmonic_constant + "\n")

