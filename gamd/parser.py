'''
Created on Oct 28, 2020

Define the Config, Parser, and ParserFactory objects for reading and storing
GaMD simulation parameters.

@author: lvotapka
'''

from __future__ import absolute_import
import xml.etree.ElementTree as ET
from xml.dom import minidom
from abc import ABCMeta, ABC
from abc import abstractmethod
from simtk import unit
from gamd import config


def strBool(bool_str):
    """
    Take the string "true" or "false" of any case and returns a 
    boolean object.
    """
    if bool_str.lower() == "true":
        return True
    elif bool_str.lower() == "false":
        return False
    else:
        raise Exception(
            "argument for strBool must be string either 'True' or 'False'.")


class Parser:
    def __init__(self):
        self.config = config.Config()
    
    def assign_value(self, value, func, useunit=None):
        if value is not None:
            if useunit is None:
                return func(value)
            else:
                return unit.Quantity(func(value), useunit)
        else:
            return None
    
    @abstractmethod
    def parse_file(self, filename):
        raise NotImplementedError("must implement parse_file")
    
    
class XmlParser(Parser):
    
    def __init__(self):
        super(XmlParser, self).__init__()
        
    def assign_tag(self, tag, func, useunit=None):
        if tag is not None:
            return self.assign_value(tag.text, func, useunit)
        else:
            return None
    
    def parse_file(self, filename):
        tree = ET.parse(filename)
        root = tree.getroot()
        for tag in root:
            if tag.tag == "temperature":
                self.config.temperature = self.assign_tag(
                    tag, float, useunit=unit.kelvin)
            elif tag.tag == "system":
                for system_tag in tag:
                    if system_tag.tag == "nonbonded-method":
                        self.config.system.nonbonded_method \
                            = self.assign_tag(system_tag, str).lower()
                    elif system_tag.tag == "nonbonded-cutoff":
                        self.config.system.nonbonded_cutoff \
                            = self.assign_tag(system_tag, float, 
                                              useunit=unit.nanometer)
                    elif system_tag.tag == "constraints":
                        self.config.system.constraints \
                            = self.assign_tag(system_tag, str).lower()
                    else:
                        print("Warning: parameter in XML not found in system "\
                              "tag. Spelling error?", system_tag.tag)
            elif tag.tag == "barostat":
                self.config.barostat = config.BarostatConfig()
                for barostat_tag in tag:
                    if barostat_tag.tag == "pressure":
                        self.config.barostat.pressure \
                            = self.assign_tag(barostat_tag, float, 
                                              useunit=unit.bar)
                    elif barostat_tag.tag == "frequency":
                        self.config.barostat.frequency \
                            = self.assign_tag(barostat_tag, int)
                    else:
                        print("Warning: parameter in XML not found in "\
                              "barostat tag. Spelling error?", barostat_tag.tag)
            elif tag.tag == "run-minimization":
                self.config.run_minimization = self.assign_tag(
                    tag, strBool)
            elif tag.tag == "integrator":
                for integrator_tag in tag:
                    if integrator_tag.tag == "algorithm":
                        self.config.integrator.algorithm \
                            = self.assign_tag(integrator_tag, str).lower()
                    elif integrator_tag.tag == "boost-type":
                        self.config.integrator.boost_type \
                            = self.assign_tag(integrator_tag, str).lower()
                    elif integrator_tag.tag == "sigma0":
                        for sigma0_tag in integrator_tag:
                            if sigma0_tag.tag == "primary":
                                self.config.integrator.sigma0.primary \
                                    = self.assign_tag(
                                        sigma0_tag, float, 
                                        useunit=unit.kilocalories_per_mole)
                            elif sigma0_tag.tag == "secondary":
                                self.config.integrator.sigma0.secondary \
                                    = self.assign_tag(
                                        sigma0_tag, float, 
                                        useunit=unit.kilocalories_per_mole)
                            else:
                                print("Warning: parameter in XML not found in "\
                                      "sigma0 tag. Spelling error?", 
                                      sigma0_tag.tag)
                    elif integrator_tag.tag == "random-seed":
                        self.config.integrator.random_seed \
                            = self.assign_tag(integrator_tag, int)
                    elif integrator_tag.tag == "dt":
                        self.config.integrator.dt \
                            = self.assign_tag(integrator_tag, float, 
                                              useunit=unit.picoseconds)
                    elif integrator_tag.tag == "friction-coefficient":
                        self.config.integrator.friction_coefficient \
                            = self.assign_tag(integrator_tag, float, 
                                              useunit=unit.picoseconds**-1)
                    elif integrator_tag.tag == "number-of-steps":
                        for number_steps_tag in integrator_tag:
                            if number_steps_tag.tag == "conventional-md-prep":
                                self.config.integrator.number_of_steps\
                                    .conventional_md_prep \
                                    = self.assign_tag(number_steps_tag, int)
                            elif number_steps_tag.tag == "conventional-md":
                                self.config.integrator.number_of_steps\
                                    .conventional_md \
                                    = self.assign_tag(number_steps_tag, int)
                            elif number_steps_tag.tag == "gamd-equilibration-prep":
                                self.config.integrator.number_of_steps\
                                    .gamd_equilibration_prep \
                                    = self.assign_tag(number_steps_tag, int)
                            elif number_steps_tag.tag == "gamd-equilibration":
                                self.config.integrator.number_of_steps\
                                    .gamd_equilibration \
                                    = self.assign_tag(number_steps_tag, int)
                            elif number_steps_tag.tag == "gamd-production":
                                self.config.integrator.number_of_steps\
                                    .gamd_production \
                                    = self.assign_tag(number_steps_tag, int)
                            elif number_steps_tag.tag == "total-simulation-length":
                                self.config.integrator.number_of_steps\
                                    .total_simulation_length \
                                    = self.assign_tag(number_steps_tag, int)
                            elif number_steps_tag.tag == "averaging-window-interval":
                                self.config.integrator.number_of_steps\
                                    .averaging_window_interval \
                                    = self.assign_tag(number_steps_tag, int)
                            else:
                                print("Warning: parameter in XML not found in "\
                                      "number-of-steps tag. Spelling error?", 
                                      number_steps_tag.tag)
                    else:
                        print("Warning: parameter in XML not found in "\
                              "integrator tag. Spelling error?", 
                              integrator_tag.tag)
            elif tag.tag == "input-files":
                for input_files_tag in tag:
                    if input_files_tag.tag == "amber":
                        amber_config = config.AmberConfig()
                        for amber_tag in input_files_tag:
                            if amber_tag.tag == "topology":
                                amber_config.topology = self.assign_tag(
                                    amber_tag, str)
                            elif amber_tag.tag == "coordinates":
                                amber_config.coordinates = self.assign_tag(
                                    amber_tag, str)
                                if "type" in amber_tag.attrib:
                                    type_attrib = amber_tag.attrib["type"]
                                    amber_config.coordinates_filetype \
                                        = type_attrib
                            else:
                                print("Warning: parameter in XML not found in "\
                                      "amber tag. Spelling error?", 
                                      amber_tag.tag)
                        self.config.input_files.amber = amber_config
                        
                    elif input_files_tag.tag == "charmm":
                        charmm_config = config.CharmmConfig()
                        for charmm_tag in input_files_tag:
                            if charmm_tag.tag == "topology":
                                charmm_config.topology = self.assign_tag(
                                    charmm_tag, str)
                            elif charmm_tag.tag == "coordinates":
                                charmm_config.coordinates = self.assign_tag(
                                    charmm_tag, str)
                            elif charmm_tag.tag == "parameters":
                                charmm_config.parameters = []
                                for xml_params_filename in charmm_tag:
                                    charmm_config.parameters.append(
                                        self.assign_tag(xml_params_filename, 
                                                        str))
                            else:
                                print("Warning: parameter in XML not found in "\
                                      "charmm tag. Spelling error?", 
                                      charmm_tag.tag)
                        self.config.input_files.charmm = charmm_config
                    elif input_files_tag.tag == "gromacs":
                        gromacs_config = config.GromacsConfig()
                        for gromacs_tag in input_files_tag:
                            if gromacs_tag.tag == "topology":
                                gromacs_config.topology = self.assign_tag(
                                    gromacs_tag, str)
                            elif gromacs_tag.tag == "coordinates":
                                gromacs_config.coordinates = \
                                    self.assign_tag(gromacs_tag, str)
                            elif gromacs_tag.tag == "include-dir":
                                gromacs_config.include_dir = \
                                    self.assign_tag(gromacs_tag, str)
                            else:
                                print("Warning: parameter in XML not "\
                                      "found in gromacs tag. Spelling "\
                                      "error?", gromacs_tag.tag)
                        self.config.input_files.gromacs = gromacs_config
                        
                    elif input_files_tag.tag == "forcefield":
                        forcefield_config = config.ForceFieldConfig()
                        for forcefield_tag in input_files_tag:
                            if forcefield_tag.tag == "coordinates":
                                forcefield_config.coordinates = self.assign_tag(
                                    forcefield_tag, str)
                                    
                            elif forcefield_tag.tag == "forcefields":
                                for forcefields_tag in forcefield_tag:
                                    if forcefields_tag.tag == "native":
                                        forcefield_config\
                                            .forcefield_list_native = []
                                        for file_tag in forcefields_tag:
                                            forcefield_config\
                                                .forcefield_list_native\
                                                .append(self.assign_tag(
                                                    file_tag, str))
                                    elif forcefields_tag.tag == "external":
                                        forcefield_config\
                                            .forcefield_list_external = []
                                        for file_tag in forcefields_tag:
                                            forcefield_config\
                                            .forcefield_list_external\
                                            .append(self.assign_tag(
                                                file_tag, str))
                                    else:
                                        print("Warning: parameter in XML not "\
                                              "found in forcefields tag. "\
                                              "Spelling error?", 
                                              forcefields_tag.tag)
                                    
                            else:
                                print("Warning: parameter in XML not found in "\
                                      "forcefield tag. Spelling error?", 
                                      forcefield_tag.tag)
                        
                        self.config.input_files.forcefield = forcefield_config
                
                    else:
                        raise Exception("input-files type not implemented:", 
                                        input_files_tag.tag)
            
            elif tag.tag == "outputs":
                for outputs_tag in tag:
                    if outputs_tag.tag == "directory":
                        self.config.outputs.directory \
                            = self.assign_tag(outputs_tag, str)
                    elif outputs_tag.tag == "overwrite-output":
                        self.config.outputs.overwrite_output \
                            = self.assign_tag(outputs_tag, strBool)
                    elif outputs_tag.tag == "reporting":
                        for reporting_tag in outputs_tag:
                            if reporting_tag.tag == "energy":
                                for energy_tag in reporting_tag:
                                    if energy_tag.tag == "interval":
                                        self.config.outputs.reporting\
                                            .energy_interval \
                                            = self.assign_tag(energy_tag, int)
                                    else:
                                        print("Warning: parameter in XML not "\
                                            "found in energy tag. Spelling "\
                                            "error?", energy_tag.tag)
                            elif reporting_tag.tag == "coordinates":
                                for coordinates_tag in reporting_tag:
                                    if coordinates_tag.tag == "interval":
                                        self.config.outputs.reporting\
                                            .coordinates_interval \
                                            = self.assign_tag(coordinates_tag, 
                                                              int)
                                    elif coordinates_tag.tag == "file-type":
                                        self.config.outputs.reporting\
                                            .coordinates_file_type \
                                            = self.assign_tag(coordinates_tag, 
                                                              str).lower()
                                    else:
                                        print("Warning: parameter in XML not "\
                                            "found in coordinates tag. "\
                                            "Spelling error?", 
                                            coordinates_tag.tag)
                            elif reporting_tag.tag == "restart-checkpoint":
                                for restart_tag in reporting_tag:
                                    if restart_tag.tag == "interval":
                                        self.config.outputs.reporting\
                                            .restart_checkpoint_interval \
                                            = self.assign_tag(restart_tag, int)
                                    else:
                                        print("Warning: parameter in XML not "\
                                            "found in restart-checkpoint tag. "\
                                            "Spelling error?", restart_tag.tag)
                            elif reporting_tag.tag == "statistics":
                                for statistics_tag in reporting_tag:
                                    if statistics_tag.tag == "interval":
                                        self.config.outputs.reporting\
                                            .statistics_interval \
                                            = self.assign_tag(statistics_tag, 
                                                              int)
                                    else:
                                        print("Warning: parameter in XML not "\
                                            "found in statistics tag. "\
                                            "Spelling error?", 
                                            statistics_tag.tag)
                            else:
                                print("Warning: parameter in XML not found in "\
                                      "reporting tag. Spelling error?", 
                                      reporting_tag.tag)
                    else:
                        print("Warning: parameter in XML not found in "\
                              "outputs tag. Spelling error?", 
                              outputs_tag.tag)
            
            else:
                print("Warning: parameter in XML not found in config. "\
                      "Spelling error?", tag.tag)


class ParserFactory:
    def __init__(self):
        return
        
    def parse_file(self, input_file, input_type):
        input_type = input_type.lower()
        if input_type == "xml":
            myparser = XmlParser()
            myparser.parse_file(input_file)
            config = myparser.config
        else:
            raise Exception("input type not implemented: %s" % input_type)
        
        return config


if __name__ == "__main__":
    myparser = XmlParser()
    myparser.parse_file("../example.xml")
    myparser.config.serialize("/tmp/gamdconfig.xml")