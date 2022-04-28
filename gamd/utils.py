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
    def __init__(self, filename, mode, integrator):
        self.gamdDatFile = None
        self.filename = filename
        self.tracked_integrator = integrator
        self.tracked_names = ["Vmax", "Vmin", "Vavg", "sigmaV",
                              "k0", "k", "sigma0", "threshold_energy"]
        self.tracked_values = self.TrackedValues(integrator, self.tracked_names)
        self.headers = self.tracked_values.get_names()

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

        for header in self.headers:
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
        if self.tracked_values.did_values_change():
            step = simulation.currentStep
            self.tracked_values.update_values()
            output_string = self.__create_output_row(step)
            self.gamdDatFile.write(output_string + "\n")



    def __create_output_row(self, step):
        output_string = str(step)
        values = self.tracked_values.get_values()
        for header in self.headers:
            output_string = output_string + ", " + str(values[header])
        return output_string

    class TrackedValues:
        def __init__(self, tracked_integrator, names):
            self.__integrator = tracked_integrator
            self.names = names
            self.tracked_names = []
            for name in self.names:
                actual_names = tracked_integrator.get_names(name)
                for actual_name in actual_names:
                    self.tracked_names.append(actual_name)

            self.tracked_values = {}
            self.update_values()

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
