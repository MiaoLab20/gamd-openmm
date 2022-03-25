"""
gamd.py: Implements the GaMD integration method.

Portions copyright (c) 2020 University of Kansas
Authors: Matthew Copeland

"""

from openmm.app.statedatareporter import StateDataReporter
import openmm.unit as unit
from .stage_integrator import BoostType


def create_gamd_log(gamdLog, filename):
    with open(filename, 'w') as f:
        keys = list(gamdLog[0])
        for header in keys[:-1]:
            f.write(header + ", ")
        f.write(keys[-1] + "\n")
        for entry in gamdLog:
            for header in keys[:-1]:
                f.write(str(entry[header]) + ", ")
            f.write(str(entry[keys[-1]]) + "\n")


class ExpandedStateDataReporter(StateDataReporter):

    def __init__(self, system, file, reportInterval, step=False,
                 time=False, brokenOutForceEnergies=False,
                 potentialEnergy=False, kineticEnergy=False, totalEnergy=False,
                 temperature=False, volume=False, density=False,
                 progress=False, remainingTime=False, speed=False,
                 elapsedTime=False, separator=',', systemMass=None,
                 totalSteps=None):

        self._brokenOutForceEnergies = brokenOutForceEnergies
        self._system = system
        super().__init__(file, reportInterval, step, time,
                         potentialEnergy, kineticEnergy, totalEnergy,
                         temperature, volume, density, progress, remainingTime,
                         speed, elapsedTime, separator, systemMass, totalSteps)

    def _constructReportValues(self, simulation, state):
        values = super()._constructReportValues(simulation, state)
        if self._brokenOutForceEnergies:
            for i, force in enumerate(self._system.getForces()):
                values.append(simulation.context.getState(
                    getEnergy=True,
                    groups={i}).getPotentialEnergy().value_in_unit(
                    unit.kilojoules_per_mole))
        return values

    def _constructHeaders(self):
        headers = super()._constructHeaders()
        if self._brokenOutForceEnergies:
            for i, force in enumerate(self._system.getForces()):
                headers.append(force.__class__.__name__)

        return headers


class GamdDatReporter:
    def __init__(self, filename, mode, integrator,
                 first_boost_type, first_boost_group,
                 second_boost_type, second_boost_group):
        self.gamdDatFile = None
        self.filename = filename
        self.tracked_integrator = integrator
        self.tracked_values = []
        if (first_boost_type == BoostType.DUAL_TOTAL_DIHEDRAL or
                second_boost_type == BoostType.DUAL_TOTAL_DIHEDRAL or
                first_boost_type == BoostType.DUAL_NON_BONDED_DIHEDRAL or
                second_boost_type == BoostType.DUAL_NON_BONDED_DIHEDRAL):
            raise ValueError("The GamdLogger expects single value boost types as arguments, not compound boost types."
                             " Compound boost types should be broken up.\n")

        if second_boost_type == BoostType.TOTAL:
            self.tracked_values.append(self.TrackedValue(second_boost_type, second_boost_group, integrator))
            self.tracked_values.append(self.TrackedValue(first_boost_type, first_boost_group, integrator))
        else:
            self.tracked_values.append(self.TrackedValue(first_boost_type, first_boost_group, integrator))
            self.tracked_values.append(self.TrackedValue(second_boost_type, second_boost_group, integrator))

        # write_headers uses tracked_values, so we want to only write it out after tracked values is set.
        self.gamdDatFile = open(filename, mode)
        if mode == "w":
            self.__write_header()
            output_string = self.__create_output_row(1)
            self.gamdDatFile.write(output_string + "\n")

    def __del__(self):
        if self.gamdDatFile is not None:
            self.close()

    def close(self):
        self.gamdDatFile.close()
        self.gamdDatFile = None

    def __write_header(self):
        header_str = "step"
        for tracked_value in self.tracked_values:
            headers = tracked_value.get_names()
            for header in headers:
                header_str = header_str + ", " + header
        self.gamdDatFile.write(header_str + "\n")

    def describeNextReport(self, simulation):
        """
            We are basically fudging this one, since we are using the reporter to track the changes
            and print out an update, when a tracked value has changed.
            The returned tuple is the number of steps until we want to be executed again, and then
            whether we want positions, velocities, forces or energies reported to us.

            Unfortunately, we don't get the GaMD variables this way.
        """
        steps = 1
        return steps, False, False, False, True

    def report(self, simulation, state):
        changes_occurred = False
        for tracked_value_set in self.tracked_values:
            changes_occurred = changes_occurred or tracked_value_set.did_values_change()
        if changes_occurred:
            step = simulation.currentStep
            output_string = self.__create_output_row(step)
            self.gamdDatFile.write(output_string + "\n")

            for tracked_value_set in self.tracked_values:
                tracked_value_set.update_values()

    def __create_output_row(self, step):
        output_string = str(step)
        for tracked_value_set in self.tracked_values:
            # We want to make sure that our values get printed out in the same order
            # as our headers did, rather than whatever order a "for key, value" would
            # give us.
            values = tracked_value_set.get_values()
            headers = tracked_value_set.get_names()
            for header in headers:
                output_string = output_string + ", " + str(values[header])
        return output_string

    class TrackedValue:
        def __init__(self, boost_type, group, tracked_integrator):
            self.__boost_type = boost_type
            self.__group = group
            self.__integrator = tracked_integrator
            self.tracked_names = [tracked_integrator.get_variable_name_by_type(boost_type, "Vmax"),
                                  tracked_integrator.get_variable_name_by_type(boost_type, "Vmin"),
                                  tracked_integrator.get_variable_name_by_type(boost_type, "Vavg"),
                                  tracked_integrator.get_variable_name_by_type(boost_type, "sigmaV"),
                                  tracked_integrator.get_variable_name_by_type(boost_type, "k0"),
                                  tracked_integrator.get_variable_name_by_type(boost_type, "k"),
                                  tracked_integrator.get_variable_name_by_type(boost_type, "sigma0"),
                                  tracked_integrator.get_variable_name_by_type(boost_type, "threshold_energy")]
            self.tracked_values = {}
            self.update_values()

        def get_boost_type(self):
            return self.__boost_type

        def did_values_change(self):
            result = False
            for name in self.tracked_names:
                if self.__integrator.getGlobalVariableByName(name) != self.tracked_values[name]:
                    result = True
                    break
            return result

        def update_values(self):
            for name in self.tracked_names:
                self.tracked_values[name] = self.__integrator.getGlobalVariableByName(name)

        def get_names(self):
            return self.tracked_names

        def get_values(self):
            return self.tracked_values
