
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import sys

class DebugLogger:

    def __init__(self, filename, mode):
        self.filename = filename
        self.debugLog = open(filename, mode)

    def __del__(self):
        self.debugLog.close()

    def write_global_variables_headers(self, integrator):
        number_of_globals = integrator.getNumGlobalVariables()
        second_to_last = number_of_globals - 1
        for index in range(0, number_of_globals):
            name = integrator.getGlobalVariableName(index)
            self.debugLog.write(str(name))
            if index < second_to_last:
                self.debugLog.write(",")
        self.debugLog.write("\n")

    def write_global_variables_values(self, integrator):
        number_of_globals = integrator.getNumGlobalVariables()
        second_to_last = number_of_globals - 1
        for index in range(0, number_of_globals):
            name = integrator.getGlobalVariableName(index)
            value = integrator.getGlobalVariableByName(name)
            self.debugLog.write(str(value))
            if index < second_to_last:
                self.debugLog.write(",")
        self.debugLog.write("\n")

    @staticmethod
    def print_integration_algorithm_to_screen(integrator):
        for i in range(integrator.getNumComputations()):
            print(integrator.getComputationStep(i))
        sys.exit(-1)

    @staticmethod
    def print_global_variables_to_screen(integrator):
        for index in range(0, integrator.getNumGlobalVariables()):
            name = integrator.getGlobalVariableName(index)
            value = integrator.getGlobalVariableByName(name)
            print(name + ":  " + str(value))
