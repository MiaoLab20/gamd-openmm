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
            if tag.tag == "system_files":
                xml_system_files_type_text = tag.find("type").text.lower()
                if xml_system_files_type_text == "amber":
                    amber_config = config.AmberConfig()
                    for amber_tag in tag:
                        if amber_tag.tag == "prmtop_filename":
                            amber_config.prmtop_filename = \
                                self.assign_tag(amber_tag, str)
                                
                        elif amber_tag.tag == "inpcrd_filename":
                            amber_config.inpcrd_filename = \
                                self.assign_tag(amber_tag, str)
                                
                        elif amber_tag.tag == \
                                "load_box_vectors_from_coordinates_file":
                            amber_config.load_box_vecs_from_coords_file = \
                                self.assign_tag(amber_tag, strBool)
                        
                        elif amber_tag.tag == "type":
                            pass
                        
                        else:
                            print("Warning: parameter in XML not found in "\
                                  "amber_config. Spelling error?", 
                                  amber_tag.tag)
                    
                    self.config.system_files_config = amber_config
                
                elif xml_system_files_type_text == "charmm":
                    charmm_config = config.CharmmConfig()
                    for charmm_tag in tag:
                        if charmm_tag.tag == "psf_filename":
                            charmm_config.psf_filename = \
                                self.assign_tag(charmm_tag, str)
                                
                        elif charmm_tag.tag == "pdb_filename":
                            charmm_config.pdb_filename = \
                                self.assign_tag(charmm_tag, str)
                                
                        elif charmm_tag.tag == \
                                "params_filenames":
                            charmm_config.params_filenames = []
                            for xml_params_filename in charmm_tag:
                                charmm_config.params_filenames.append(
                                    self.assign_tag(xml_params_filename, str))
                        
                        elif charmm_tag.tag == "type":
                            pass
                        
                        else:
                            print("Warning: parameter in XML not found in "\
                                  "charmm_config. Spelling error?", 
                                  charmm_tag.tag)
                    
                    self.config.system_files_config = charmm_config
                
                elif xml_system_files_type_text == "gromacs":
                    gromacs_config = config.GromacsConfig()
                    for gro_tag in tag:
                        if gro_tag.tag == "gro_filename":
                            gromacs_config.gro_filename = \
                                self.assign_tag(gro_tag, str)
                                
                        elif gro_tag.tag == "top_filename":
                            gromacs_config.top_filename = \
                                self.assign_tag(gro_tag, str)
                                
                        elif gro_tag.tag == "include_dir":
                            gromacs_config.include_dir = \
                                self.assign_tag(gro_tag, str)
                        
                        elif gro_tag.tag == "type":
                            pass
                        
                        else:
                            print("Warning: parameter in XML not found in "\
                                  "gromacs_config. Spelling error?", 
                                  gro_tag.tag)
                    
                    self.config.system_files_config = gromacs_config
                
                elif xml_system_files_type_text == "forcefield":
                    forcefield_config = config.ForceFieldConfig()
                    for forcefield_tag in tag:
                        if forcefield_tag.tag == "pdb_filename":
                            forcefield_config.pdb_filename = \
                                self.assign_tag(forcefield_tag, str)
                                
                        elif forcefield_tag.tag == \
                                "forcefield_filenames":
                            forcefield_config.forcefield_filenames = []
                            for forcefield_filename in forcefield_tag:
                                forcefield_config.forcefield_filenames.append(
                                    self.assign_tag(forcefield_filename, str))
                        
                        elif forcefield_config.tag == "type":
                            pass
                        
                        else:
                            print("Warning: parameter in XML not found in "\
                                  "charmm_config. Spelling error?", 
                                  forcefield_config.tag)
                    
                    self.config.system_files_config = forcefield_config
                
                else:
                    raise Exception("system_files type not implemented:", 
                                    xml_system_files_type_text)
            
            elif tag.tag == "box_vectors":
                config.box_vectors = config.deserialize_box_vectors(tag)
            
            elif tag.tag == "output_directory":
                self.config.output_directory = self.assign_tag(tag, str)
                
            elif tag.tag == "overwrite_output":
                self.config.overwrite_output = self.assign_tag(tag, strBool)
                
            elif tag.tag == "chunk_size":
                self.config.chunk_size = self.assign_tag(tag, int)
            
            elif tag.tag == "nonbonded_method":
                self.config.nonbonded_method = \
                    self.assign_tag(tag, str).lower()
                
            elif tag.tag == "nonbonded_cutoff":
                self.config.nonbonded_cutoff = self.assign_tag(
                    tag, float, useunit=unit.nanometer)
                
            elif tag.tag == "constraints":
                self.config.constraints = self.assign_tag(tag, str).lower()
            
            elif tag.tag == "integrator_type":
                self.config.integrator_type = self.assign_tag(tag, str).lower()
            
            elif tag.tag == "friction_coefficient":
                self.config.friction_coefficient = self.assign_tag(
                    tag, float, useunit=unit.picosecond**-1)
                
            elif tag.tag == "target_temperature":
                self.config.target_temperature = self.assign_tag(
                    tag, float, useunit=unit.kelvin)
                
            elif tag.tag == "random_seed":
                self.config.random_seed = self.assign_tag(tag, int)
                
            elif tag.tag == "dt":
                self.config.dt = self.assign_tag(
                    tag, float, useunit=unit.picosecond)
                
            elif tag.tag == "use_barostat":
                self.config.use_barostat = self.assign_tag(tag, strBool)
            
            elif tag.tag == "barostat_target_pressure":
                self.config.barostat_target_pressure = self.assign_tag(
                    tag, float, useunit=unit.bar)
                
            elif tag.tag == "barostat_target_temperature":
                self.config.barostat_target_temperature = self.assign_tag(
                    tag, float, useunit=unit.kelvin)
            
            elif tag.tag == "barostat_frequency":
                self.config.barostat_frequency = self.assign_tag(tag, int)
            
            elif tag.tag == "run_minimization":
                self.config.run_minimization = self.assign_tag(tag, strBool)
                
            elif tag.tag == "initial_temperature":
                self.config.initial_temperature = self.assign_tag(
                    tag, float, useunit=unit.kelvin)
                
            elif tag.tag == "energy_reporter_frequency":
                self.config.energy_reporter_frequency = self.assign_tag(
                    tag, int)
                
            elif tag.tag == "coordinates_reporter_frequency":
                self.config.coordinates_reporter_frequency = self.assign_tag(
                    tag, int)
                
            elif tag.tag == "coordinates_reporter_file_type":
                self.config.coordinates_reporter_file_type = self.assign_tag(
                    tag, str).lower()
            
            elif tag.tag == "gamd_bound":
                self.config.gamd_bound = self.assign_tag(tag, str).lower()
            
            elif tag.tag == "total_simulation_length":
                self.config.total_simulation_length = self.assign_tag(tag, int)
            
            elif tag.tag == "total_boost":
                self.config.total_boost = self.assign_tag(tag, strBool)
                
            elif tag.tag == "total_boost_sigma0":
                self.config.total_boost_sigma0 = self.assign_tag(
                    tag, float, useunit=unit.kilocalories_per_mole)
            
            elif tag.tag == "dihedral_boost":
                self.config.dihedral_boost = self.assign_tag(tag, strBool)
                
            elif tag.tag == "dihedral_boost_sigma0":
                self.config.dihedral_boost_sigma0 = self.assign_tag(
                    tag, float, useunit=unit.kilocalories_per_mole)
            
            elif tag.tag == "num_steps_conventional_md":
                self.config.num_steps_conventional_md = self.assign_tag(
                    tag, int)
                
            elif tag.tag == "num_steps_conventional_md_prep":
                self.config.num_steps_conventional_md_prep = self.assign_tag(
                    tag, int)
                
            elif tag.tag == "num_steps_per_averaging":
                self.config.num_steps_per_averaging = self.assign_tag(
                    tag, int)
                
            elif tag.tag == "num_steps_gamd_equilibration":
                self.config.num_steps_gamd_equilibration = self.assign_tag(
                    tag, int)
                
            elif tag.tag == "num_steps_gamd_equilibration_prep":
                self.config.num_steps_gamd_equilibration_prep = self.assign_tag(
                    tag, int)
            
            elif tag.tag == "restart_checkpoint_filename":
                self.config.restart_checkpoint_filename = self.assign_tag(
                    tag, str)
                
            elif tag.tag == "restart_checkpoint_frequency":
                self.config.restart_checkpoint_frequency = self.assign_tag(
                    tag, int)
            
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
            raise Exception("input type not implemented: %s", input_type)
        
        return config


if __name__ == "__main__":
    myparser = XmlParser()
    myparser.parse_file("sample_input.xml")
        


