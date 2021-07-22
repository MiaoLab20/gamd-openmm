
"""
gamd.py: Implements the GaMD integration method.

Portions copyright (c) 2020 University of Kansas
Authors: Matthew Copeland

"""

try:
    from openmm.app.statedatareporter import StateDataReporter
except:
    from simtk.openmm.app.statedatareporter import StateDataReporter
    
import simtk.unit as unit


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
        values = super()._constructReportValues(simulation,state)
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

