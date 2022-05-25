import sys
from abc import ABC
from abc import abstractmethod


class BaseDebugLogger(ABC):

    @abstractmethod
    def close(self):
        raise NotImplementedError("must implement close")

    @abstractmethod
    def write_global_variables_headers(self, integrator):
        raise NotImplementedError("must implement write_global_variables_header")

    @abstractmethod
    def write_global_variables_values(self, integrator):
        raise NotImplementedError("must implement write_global_variables_values")

    @staticmethod
    def print_integration_algorithm_to_screen(integrator):
        for i in range(integrator.getNumComputations()):
            print(integrator.getComputationStep(i))
        sys.exit(-1)

    @staticmethod
    def write_integration_algorithm_to_file(filename, integrator):
        with open(filename, "a") as integration_algo_file:
            integration_algo_file.write(integrator.__class__.__module__ + "." + integrator.__class__.__name__)
            integration_algo_file.write("\n\n")
            for i in range(integrator.getNumComputations()):
                integration_algo_file.write(str(integrator.getComputationStep(i)))
                integration_algo_file.write("\n")
            integration_algo_file.write("---------------------------------------------------------\n")

    @staticmethod
    def print_global_variables_to_screen(integrator):
        for index in range(0, integrator.getNumGlobalVariables()):
            name = integrator.getGlobalVariableName(index)
            value = integrator.getGlobalVariableByName(name)
            print(name + ":  " + str(value))


class NoOpDebugLogger(BaseDebugLogger):
    """
    The intent of this class is to provide an object that
    will act like a DebugLogger, but not actually write anything
    out or perform any of the expensive queries.
    """
    def __init__(self):
        pass

    def close(self):
        pass

    def write_global_variables_headers(self, integrator):
        pass

    def write_global_variables_values(self, integrator):
        pass


class DebugLogger(BaseDebugLogger):

    def __init__(self, filename, mode, denyList=None):
        self.filename = filename
        self.debugLog = open(filename, mode)
        if denyList is None:
            self.denyList = []
        else:
            self.denyList = denyList

    def __del__(self):
        self.debugLog.close()

    def close(self):
        self.debugLog.close()

    def __get_all_headers(self, integrator):
        result = []
        number_of_globals = integrator.getNumGlobalVariables()
        for index in range(0, number_of_globals):
            name = integrator.getGlobalVariableName(index)
            result.append(name)
        return result

    def __get_filtered_headers(self, integrator):
        all_headers = self.__get_all_headers(integrator)
        headers = [header for header in all_headers if header not in self.denyList]
        return headers

    def write_global_variables_headers(self, integrator):
        headers = self.__get_filtered_headers(integrator)
        headers_string = ",".join(map(str, headers))
        self.debugLog.write(headers_string)
        self.debugLog.write("\n")

    def write_global_variables_values(self, integrator):
        headers = self.__get_filtered_headers(integrator)
        values = []
        for header in headers:
            values.append(integrator.getGlobalVariableByName(header))
        values_string = ",".join(map(str, values))
        self.debugLog.write(str(values_string))
        self.debugLog.write("\n")

