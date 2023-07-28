"""
Created on Oct 29, 2020

Represents a configuration for a GaMD run: all the settings provided by the
user.

@author: lvotapka
"""

import xml.etree.ElementTree as ET
from xml.dom import minidom

import numpy as np
import openmm.unit as unit


def assign_tag(root, tagname, value, attributes=None):
    if attributes is None:
        attributes = {}
    xmlTag = ET.SubElement(root, tagname)
    if value is not None:
        xmlTag.text = str(value)
    for attribute in attributes:
        xmlTag.set(attribute, attributes[attribute])

    return


class SystemConfig:
    def __init__(self):
        self.nonbonded_method = "PME"
        self.nonbonded_cutoff = 0.9*unit.nanometer
        self.constraints = "HBonds"
        self.switch_distance = 1.0*unit.nanometer
        self.ewald_error_tolerance = 0.0005
        return

    def serialize(self, root):
        assign_tag(root, "nonbonded-method", self.nonbonded_method)
        assign_tag(root, "nonbonded-cutoff", self.nonbonded_cutoff.value_in_unit(unit.nanometers))
        assign_tag(root, "constraints", self.constraints)
        assign_tag(root, "switch-distance", self.switch_distance.value_in_unit(unit.nanometers))
        assign_tag(root, "ewald-error-tolerance", self.ewald_error_tolerance)
        return


class BarostatConfig:
    def __init__(self):
        self.pressure = 1.0 * unit.bar
        self.frequency = 25
        return

    def serialize(self, root):
        assign_tag(root, "pressure", self.pressure.value_in_unit(unit.bar))
        assign_tag(root, "frequency", self.frequency)
        return


class IntegratorSigmaConfig:
    def __init__(self):
        self.primary = 6.0 * unit.kilocalories_per_mole
        self.secondary = 6.0 * unit.kilocalories_per_mole
        return

    def serialize(self, root):
        assign_tag(root, "primary", self.primary.value_in_unit(unit.kilocalories_per_mole))
        assign_tag(root, "secondary", self.secondary.value_in_unit(unit.kilocalories_per_mole))
        return


class IntegratorNumberOfStepsConfig:
    def __init__(self):
        self.conventional_md_prep = 0
        self.conventional_md = 0
        self.gamd_equilibration_prep = 0
        self.gamd_equilibration = 0
        self.gamd_production = 0
        self.total_simulation_length = 0
        self.averaging_window_interval = 0
        return

    def serialize(self, root):
        assign_tag(root, "conventional-md-prep", self.conventional_md_prep)
        assign_tag(root, "conventional-md", self.conventional_md)
        assign_tag(root, "gamd-equilibration-prep", self.gamd_equilibration_prep)
        assign_tag(root, "gamd-equilibration", self.gamd_equilibration)
        assign_tag(root, "gamd-production", self.gamd_production)
        #assign_tag(root, "total-simulation-length", self.total_simulation_length)
        assign_tag(root, "averaging-window-interval", self.averaging_window_interval)
        return

    def compute_total_simulation_length(self):
        self.total_simulation_length = self.conventional_md \
            + self.gamd_equilibration + self.gamd_production


class IntegratorConfig:
    def __init__(self):
        self.algorithm = "langevin"
        self.boost_type = "lower-dual"
        self.sigma0 = IntegratorSigmaConfig()
        self.random_seed = 0
        self.dt = 0.002 * unit.picoseconds
        self.friction_coefficient = 1.0 * unit.picoseconds ** -1
        self.number_of_steps = IntegratorNumberOfStepsConfig()
        return

    def serialize(self, root):
        assign_tag(root, "algorithm", self.algorithm)
        assign_tag(root, "boost-type", self.boost_type)
        xml_sigma0_tags = ET.SubElement(root, "sigma0")
        self.sigma0.serialize(xml_sigma0_tags)
        assign_tag(root, "random-seed", self.random_seed)
        assign_tag(root, "dt", self.dt.value_in_unit(unit.picoseconds))
        assign_tag(root, "friction-coefficient", self.friction_coefficient.value_in_unit(unit.picoseconds**-1))
        xml_number_of_steps_tags = ET.SubElement(root, "number-of-steps")
        self.number_of_steps.serialize(xml_number_of_steps_tags)
        return


class AmberConfig:
    def __init__(self):
        self.topology = ""
        self.coordinates = ""
        self.coordinates_filetype = ""
        return

    def serialize(self, root):
        assign_tag(root, "topology", self.topology)
        assign_tag(root, "coordinates", self.coordinates, {"type": self.coordinates_filetype})
        return


class CharmmConfig:
    def __init__(self):
        self.topology = ""
        self.coordinates = ""
        self.coordinates_filetype = ""
        self.parameters = []
        self.box_vectors = []
        self.is_config_box_vector_defined = False
        return

    def serialize(self, root):
        assign_tag(root, "topology", self.topology)
        assign_tag(root, "coordinates", self.coordinates, {"type": self.coordinates_filetype})
        parameters_tag = ET.SubElement(root, "parameters")
        for parameter in self.parameters:
            assign_tag(parameters_tag, "parameters", parameter)
        box_vector_tag = ET.SubElement(root, "box-vectors")
        assign_tag(box_vector_tag, "a", self.box_vectors[0].value_in_unit(unit.nanometer))
        assign_tag(box_vector_tag, "b", self.box_vectors[1].value_in_unit(unit.nanometer))
        assign_tag(box_vector_tag, "c", self.box_vectors[2].value_in_unit(unit.nanometer))
        assign_tag(box_vector_tag, "alpha", self.box_vectors[3].value_in_unit(unit.degree))
        assign_tag(box_vector_tag, "beta", self.box_vectors[4].value_in_unit(unit.degree))
        assign_tag(box_vector_tag, "gamma", self.box_vectors[5].value_in_unit(unit.degree))

        return


class GromacsConfig:
    def __init__(self):
        self.topology = ""
        self.coordinates = ""
        self.include_dir = ""
        return

    def serialize(self, root):
        assign_tag(root, "topology", self.topology)
        assign_tag(root, "coordinates", self.coordinates)
        assign_tag(root, "include-dir", self.include_dir)
        return


class ForceFieldConfig:
    def __init__(self):
        self.coordinates = ""
        self.forcefield_list_native = []
        self.forcefield_list_external = []
        return

    def serialize(self, root):
        assign_tag(root, "coordinates", self.coordinates)
        xmlForcefields = ET.SubElement(root, "forcefields")
        xmlNative = ET.SubElement(xmlForcefields, "native")
        xmlExternal = ET.SubElement(xmlForcefields, "external")
        for native_filename in self.forcefield_list_native:
            assign_tag(xmlNative, "file", native_filename)
        for external_filename in self.forcefield_list_external:
            assign_tag(xmlExternal, "file", external_filename)
        return


class InputFilesConfig:
    def __init__(self):
        self.amber = None
        self.charmm = None
        self.gromacs = None
        self.forcefield = None
        return

    def serialize(self, root):
        if self.amber is not None:
            xml_amber_tags = ET.SubElement(root, "amber")
            self.amber.serialize(xml_amber_tags)
        if self.charmm is not None:
            xml_charmm_tags = ET.SubElement(root, "charmm")
            self.charmm.serialize(xml_charmm_tags)
        if self.gromacs is not None:
            xml_gromacs_tags = ET.SubElement(root, "gromacs")
            self.gromacs.serialize(xml_gromacs_tags)
        if self.forcefield is not None:
            xml_forcefield_tags = ET.SubElement(root, "forcefield")
            self.forcefield.serialize(xml_forcefield_tags)
        return


class OutputsReportingConfig:
    def __init__(self):
        self.energy_interval = 500
        self.coordinates_file_type = "DCD"
        self.coordinates_interval = 500
        self.restart_checkpoint_interval = 50000
        self.statistics_interval = 500
        return

    def compute_chunk_size(self):
        gcd = np.gcd.reduce(
            [self.energy_interval, self.coordinates_interval,
             self.restart_checkpoint_interval, self.statistics_interval])
        return gcd

    def serialize(self, root):
        xml_energy_tags = ET.SubElement(root, "energy")
        assign_tag(xml_energy_tags, "interval", self.energy_interval)
        xml_coordinates_tags = ET.SubElement(root, "coordinates")
        assign_tag(xml_coordinates_tags, "file-type", self.coordinates_file_type)
        assign_tag(xml_coordinates_tags, "interval", self.coordinates_interval)
        xml_statistics_tags = ET.SubElement(root, "statistics")
        assign_tag(xml_statistics_tags, "interval", self.statistics_interval)
        return


class OutputsConfig:
    def __init__(self):
        self.directory = ""
        self.overwrite_output = True
        self.reporting = OutputsReportingConfig()
        return

    def serialize(self, root):
        assign_tag(root, "directory", self.directory)
        assign_tag(root, "overwrite-output", self.overwrite_output)
        xml_reporting_tags = ET.SubElement(root, "reporting")
        self.reporting.serialize(xml_reporting_tags)
        return


class Config:
    def __init__(self):
        # set all the default values

        self.temperature = 298.15 * unit.kelvin
        self.system = SystemConfig()
        self.barostat = None #BarostatConfig()
        self.run_minimization = True
        self.integrator = IntegratorConfig()
        self.input_files = InputFilesConfig()
        self.outputs = OutputsConfig()

    def serialize(self, filename):
        root = ET.Element('gamd')
        assign_tag(root, "temperature", self.temperature.value_in_unit(unit.kelvin))
        xml_system = ET.SubElement(root, "system")
        self.system.serialize(xml_system)
        if self.barostat is not None:
            xml_barostat = ET.SubElement(root, "barostat")
            self.barostat.serialize(xml_barostat)
        assign_tag(root, "run-minimization", self.run_minimization)
        xml_integrator = ET.SubElement(root, "integrator")
        self.integrator.serialize(xml_integrator)
        xml_input_files = ET.SubElement(root, "input-files")
        self.input_files.serialize(xml_input_files)
        xml_outputs = ET.SubElement(root, "outputs")
        self.outputs.serialize(xml_outputs)

        xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(
            indent="    ")
        our_file=open(filename, 'w')
        our_file.write(xmlstr)
        our_file.close()


if __name__ == "__main__":
    pass
