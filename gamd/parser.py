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

import openmm.unit as unit

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


def assign_value(value, func, useunit=None):
    if value is not None:
        if useunit is None:
            return func(value)
        else:
            return unit.Quantity(func(value), useunit)
    else:
        return None


def assign_tag(tag, func, useunit=None):
        if tag is not None:
            return assign_value(tag.text, func, useunit)
        else:
            return None

def parse_and_assign_charmm_gui_toppar_file(charmm_config, xml_params_filename):
    extlist = ['rtf', 'prm', 'str']
    parFiles = ()
    # this input file lists a series of parameter files to be parsed. This is 
    # the default approach used by CHARMM-GUI (in their toppar.str file)
    for line in open(xml_params_filename.text, 'r'):
        if '!' in line: line = line.split('!')[0]
        parfile = line.strip()
        if len(parfile) != 0:
            ext = parfile.lower().split('.')[-1]
            if not ext in extlist: continue
            # Appending each file listed in inputfile to the existing 
            # list of parameters files to be read
            charmm_config.parameters.append(parfile)
    return charmm_config

class Parser:
    def __init__(self):
        self.config = config.Config()
    
    @abstractmethod
    def parse_file(self, filename):
        raise NotImplementedError("must implement parse_file")


def parse_system_tag(tag):
    system_config = config.SystemConfig()
    for system_tag in tag:
        if system_tag.tag == "nonbonded-method":
            system_config.nonbonded_method \
                = assign_tag(system_tag, str).lower()
        elif system_tag.tag == "nonbonded-cutoff":
            system_config.nonbonded_cutoff \
                = assign_tag(system_tag, float, 
                                  useunit=unit.nanometer)
        elif system_tag.tag == "constraints":
            system_config.constraints \
                = assign_tag(system_tag, str).lower()
        elif system_tag.tag == "switch-distance":
            system_config.switch_distance \
                = assign_tag(system_tag, float, 
                                  useunit=unit.nanometer)
        elif system_tag.tag == "ewald-error-tolerance":
            system_config.ewald_error_tolerance \
                = assign_tag(system_tag, float)
        else:
            print("Warning: parameter in XML not found in system "\
                  "tag. Spelling error?", system_tag.tag)
    return system_config


def parse_barostat_tag(tag):
    if len(tag) == 0:
        return None
        

    barostat_config = config.BarostatConfig()
    for barostat_tag in tag:
        if barostat_tag.tag == "pressure":
            barostat_config.pressure \
                = assign_tag(barostat_tag, float, 
                                  useunit=unit.bar)
        elif barostat_tag.tag == "frequency":
            barostat_config.frequency \
                = assign_tag(barostat_tag, int)
        else:
            print("Warning: parameter in XML not found in "\
                  "barostat tag. Spelling error?", barostat_tag.tag)
    return barostat_config


def parse_integrator_tag(tag):
    integrator_config = config.IntegratorConfig()
    for integrator_tag in tag:
        if integrator_tag.tag == "algorithm":
            integrator_config.algorithm  = assign_tag(
                integrator_tag, str).lower()
        elif integrator_tag.tag == "boost-type":
            integrator_config.boost_type = assign_tag(
                integrator_tag, str).lower()
        elif integrator_tag.tag == "sigma0":
            for sigma0_tag in integrator_tag:
                if sigma0_tag.tag == "primary":
                    integrator_config.sigma0.primary = assign_tag(
                        sigma0_tag, float, useunit=unit.kilocalories_per_mole)
                elif sigma0_tag.tag == "secondary":
                    integrator_config.sigma0.secondary = assign_tag(
                        sigma0_tag, float, useunit=unit.kilocalories_per_mole)
                else:
                    print("Warning: parameter in XML not found in sigma0 tag. "\
                          "Spelling error?", sigma0_tag.tag)
        elif integrator_tag.tag == "random-seed":
            integrator_config.random_seed = assign_tag(integrator_tag, int)
        elif integrator_tag.tag == "dt":
            integrator_config.dt = assign_tag(
                integrator_tag, float, useunit=unit.picoseconds)
        elif integrator_tag.tag == "friction-coefficient":
            integrator_config.friction_coefficient = assign_tag(
                integrator_tag, float, useunit=unit.picoseconds**-1)
        elif integrator_tag.tag == "number-of-steps":
            for number_steps_tag in integrator_tag:
                if number_steps_tag.tag == "conventional-md-prep":
                    integrator_config.number_of_steps.conventional_md_prep \
                        = assign_tag(number_steps_tag, int)
                elif number_steps_tag.tag == "conventional-md":
                    integrator_config.number_of_steps.conventional_md \
                        = assign_tag(number_steps_tag, int)
                elif number_steps_tag.tag == "gamd-equilibration-prep":
                    integrator_config.number_of_steps.gamd_equilibration_prep \
                        = assign_tag(number_steps_tag, int)
                elif number_steps_tag.tag == "gamd-equilibration":
                    integrator_config.number_of_steps.gamd_equilibration \
                        = assign_tag(number_steps_tag, int)
                elif number_steps_tag.tag == "gamd-production":
                    integrator_config.number_of_steps.gamd_production \
                        = assign_tag(number_steps_tag, int)
                #elif number_steps_tag.tag == "total-simulation-length":
                #    self.config.integrator.number_of_steps\
                #        .total_simulation_length \
                #        = self.assign_tag(number_steps_tag, int)
                elif number_steps_tag.tag == "averaging-window-interval":
                    integrator_config.number_of_steps\
                        .averaging_window_interval \
                        = assign_tag(number_steps_tag, int)
                else:
                    print("Warning: parameter in XML not found in "\
                          "number-of-steps tag. Spelling error?", 
                          number_steps_tag.tag)
        else:
            print("Warning: parameter in XML not found in "\
                  "integrator tag. Spelling error?", 
                  integrator_tag.tag)
            
    return integrator_config


def parse_amber_tag(input_files_tag):
    amber_config = config.AmberConfig()
    for amber_tag in input_files_tag:
        if amber_tag.tag == "topology":
            amber_config.topology = assign_tag(amber_tag, str)
        elif amber_tag.tag == "coordinates":
            amber_config.coordinates = assign_tag(amber_tag, str)
            if "type" in amber_tag.attrib:
                type_attrib = amber_tag.attrib["type"]
                amber_config.coordinates_filetype = type_attrib
        else:
            print("Warning: parameter in XML not found in amber tag. "\
                  "Spelling error?", amber_tag.tag)
    return amber_config


def parse_charmm_tag(input_files_tag):
    charmm_config = config.CharmmConfig()
    for charmm_tag in input_files_tag:
        if charmm_tag.tag == "topology":
            charmm_config.topology = assign_tag(charmm_tag, str)
        elif charmm_tag.tag == "coordinates":
            charmm_config.coordinates = assign_tag(charmm_tag, str)
            if "type" in charmm_tag.attrib:
                type_attrib = charmm_tag.attrib["type"]
                charmm_config.coordinates_filetype = type_attrib
        elif charmm_tag.tag == "box-vectors":
            box_dict = {}
            for box_info in charmm_tag:
                charmm_config.config_box_vector_defined = True
                if box_info.tag in ["a", "b", "c"]:
                    box_dict[box_info.tag] = (assign_tag(box_info, float, 
                                                        useunit=unit.nanometer))
                elif box_info.tag in ["alpha", "beta", "gamma"]:
                    box_dict[box_info.tag] = (assign_tag(box_info, float, 
                                                        useunit=unit.degree))
                else:
                    raise Exception("Unkown parameter '" +box_info.tag+ "'. "\
                         "Accepted box-vector parameters are 'a', 'b', 'c', "\
                         "'alpha', 'beta' and 'gamma'." )

            charmm_config.box_vectors = box_dict

        elif charmm_tag.tag == "parameters":
            charmm_config.parameters = []
            for xml_params_filename in charmm_tag:
                if "type" in xml_params_filename.attrib:
                    if xml_params_filename.attrib["type"] == "charmm-gui-toppar":
                        # parsing list of parameter files like CHARMM-GUI does
                        charmm_config = parse_and_assign_charmm_gui_toppar_file(
                        charmm_config, xml_params_filename)
                else:
                    charmm_config.parameters.append(
                        assign_tag(xml_params_filename, str))
        else:
            print("Warning: parameter in XML not found in "\
                  "charmm tag. Spelling error?", charmm_tag.tag)
    return charmm_config

def parse_gromacs_tag(input_files_tag):
    gromacs_config = config.GromacsConfig()
    for gromacs_tag in input_files_tag:
        if gromacs_tag.tag == "topology":
            gromacs_config.topology = assign_tag(gromacs_tag, str)
        elif gromacs_tag.tag == "coordinates":
            gromacs_config.coordinates = assign_tag(gromacs_tag, str)
        elif gromacs_tag.tag == "include-dir":
            gromacs_config.include_dir = \
                assign_tag(gromacs_tag, str)
        else:
            print("Warning: parameter in XML not found in gromacs tag. "\
                  "Spelling error?", gromacs_tag.tag)
    return gromacs_config


def parse_forcefield_tag(input_files_tag):
    forcefield_config = config.ForceFieldConfig()
    for forcefield_tag in input_files_tag:
        if forcefield_tag.tag == "coordinates":
            forcefield_config.coordinates = assign_tag(
                forcefield_tag, str)
                
        elif forcefield_tag.tag == "forcefields":
            for forcefields_tag in forcefield_tag:
                if forcefields_tag.tag == "native":
                    forcefield_config.forcefield_list_native = []
                    for file_tag in forcefields_tag:
                        forcefield_config.forcefield_list_native.append(
                            assign_tag(file_tag, str))
                elif forcefields_tag.tag == "external":
                    forcefield_config.forcefield_list_external = []
                    for file_tag in forcefields_tag:
                        forcefield_config.forcefield_list_external.append(
                            assign_tag(file_tag, str))
                else:
                    print("Warning: parameter in XML not found in forcefields "\
                          "tag. Spelling error?", forcefields_tag.tag)
                
        else:
            print("Warning: parameter in XML not found in forcefield tag. "\
                  "Spelling error?", forcefield_tag.tag)
    
    return forcefield_config


def parse_outputs_tag(tag):
    outputs_config = config.OutputsConfig()
    for outputs_tag in tag:
        if outputs_tag.tag == "directory":
            outputs_config.directory = assign_tag(outputs_tag, str)
        elif outputs_tag.tag == "overwrite-output":
            outputs_config.overwrite_output  = assign_tag(outputs_tag, strBool)
        elif outputs_tag.tag == "reporting":
            for reporting_tag in outputs_tag:
                if reporting_tag.tag == "energy":
                    for energy_tag in reporting_tag:
                        if energy_tag.tag == "interval":
                            outputs_config.reporting.energy_interval \
                                = assign_tag(energy_tag, int)
                        else:
                            print("Warning: parameter in XML not found in "\
                                  "energy tag. Spelling error?", energy_tag.tag)
                
                elif reporting_tag.tag == "coordinates":
                    for coordinates_tag in reporting_tag:
                        if coordinates_tag.tag == "file-type":
                            outputs_config.reporting.coordinates_file_type \
                                = assign_tag(coordinates_tag, str).lower()
                        else:
                            print("Warning: parameter in XML not found in "\
                                  "coordinates tag. Spelling error?", 
                                  coordinates_tag.tag)
                
                elif reporting_tag.tag == "statistics":
                    for statistics_tag in reporting_tag:
                        if statistics_tag.tag == "interval":
                            outputs_config.reporting.statistics_interval \
                                = assign_tag(statistics_tag, int)
                            outputs_config.reporting.restart_checkpoint_interval \
                                = assign_tag(statistics_tag, int)
                            outputs_config.reporting.coordinates_interval \
                                = assign_tag(statistics_tag, int)
                        
                        else:
                            print("Warning: parameter in XML not found in "\
                                  "statistics tag. Spelling error?", 
                                  statistics_tag.tag)
                else:
                    print("Warning: parameter in XML not found in "\
                          "reporting tag. Spelling error?", 
                          reporting_tag.tag)
        else:
            print("Warning: parameter in XML not found in "\
                  "outputs tag. Spelling error?", 
                  outputs_tag.tag)
            
    return outputs_config


class XmlParser(Parser):
    
    def __init__(self):
        super(XmlParser, self).__init__()
    
    def parse_file(self, filename):
        tree = ET.parse(filename)
        root = tree.getroot()
        input_file_provided = False
        for tag in root:
            if tag.tag == "temperature":
                self.config.temperature = assign_tag(
                    tag, float, useunit=unit.kelvin)
            elif tag.tag == "system":
                self.config.system = parse_system_tag(tag)
                
            elif tag.tag == "barostat":
                self.config.barostat = parse_barostat_tag(tag)
            
            elif tag.tag == "run-minimization":
                self.config.run_minimization = assign_tag(tag, strBool)
                
            elif tag.tag == "integrator":
                self.config.integrator = parse_integrator_tag(tag)
                       
            elif tag.tag == "input-files":
                for input_files_tag in tag:
                    if input_files_tag.tag == "amber":
                        assert not input_file_provided, "Only one input set "\
                            "allowed. Cannot provide more than one <amber>, "\
                            "<charmm>, <gromacs>, or <forcefield> tag."
                        
                        self.config.input_files.amber = parse_amber_tag(
                            input_files_tag)
                        input_file_provided = True
                        
                    elif input_files_tag.tag == "charmm":
                        assert not input_file_provided, "Only one input set "\
                            "allowed. Cannot provide more than one <amber>, "\
                            "<charmm>, <gromacs>, or <forcefield> tag."
                        self.config.input_files.charmm = parse_charmm_tag(
                            input_files_tag)
                        input_file_provided = True
                        
                    elif input_files_tag.tag == "gromacs":
                        assert not input_file_provided, "Only one input set "\
                            "allowed. Cannot provide more than one <amber>, "\
                            "<charmm>, <gromacs>, or <forcefield> tag."
                        self.config.input_files.gromacs = parse_gromacs_tag(
                            input_files_tag)
                        input_file_provided = True
                        
                    elif input_files_tag.tag == "forcefield":
                        assert not input_file_provided, "Only one input set "\
                            "allowed. Cannot provide more than one <amber>, "\
                            "<charmm>, <gromacs>, or <forcefield> tag."
                        self.config.input_files.forcefield = parse_forcefield_tag(
                            input_files_tag)
                        input_file_provided = True
                        
                    else:
                        raise Exception("input-files type not implemented:", 
                                        input_files_tag.tag)
            
            elif tag.tag == "outputs":
                self.config.outputs = parse_outputs_tag(tag)
            
            else:
                print("Warning: parameter in XML not found in config. "\
                      "Spelling error?", tag.tag)
        
        self.config.integrator.number_of_steps.compute_total_simulation_length()
        return


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
    pass
