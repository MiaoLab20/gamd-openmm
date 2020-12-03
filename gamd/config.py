'''
Created on Oct 29, 2020

Represents a configuration for a GaMD run: all the settings provided by the
user.

@author: lvotapka
'''

import xml.etree.ElementTree as ET
from xml.dom import minidom

from simtk import unit

def assign_tag(root, tagname, value):
    xmlTag = ET.SubElement(root, tagname)
    if value is not None:
        xmlTag.text = str(value)
    return

def serialize_box_vectors(box_vectors, xmlRoot, tagName='box_vectors'):
    '''
    Takes a 3x3 simtk.unit.Quantity representing an openMM or Parmed system box 
    vectors and serializes it as xml, saving the text into xmlRoot.
    '''
    xmlBox_vectors = ET.SubElement(xmlRoot, tagName)
    if box_vectors is not None:
        box_vectors_unitless = box_vectors.value_in_unit(unit.nanometer)
        xmlA = ET.SubElement(xmlBox_vectors, 'A')
        xmlAx = ET.SubElement(xmlA, 'x')
        xmlAx.text = str(box_vectors_unitless[0][0])
        xmlAy = ET.SubElement(xmlA, 'y')
        xmlAy.text = str(box_vectors_unitless[0][1])
        xmlAz = ET.SubElement(xmlA, 'z')
        xmlAz.text = str(box_vectors_unitless[0][2])
        xmlB = ET.SubElement(xmlBox_vectors, 'B')
        xmlBx = ET.SubElement(xmlB, 'x')
        xmlBx.text = str(box_vectors_unitless[1][0])
        xmlBy = ET.SubElement(xmlB, 'y')
        xmlBy.text = str(box_vectors_unitless[1][1])
        xmlBz = ET.SubElement(xmlB, 'z')
        xmlBz.text = str(box_vectors_unitless[1][2])
        xmlC = ET.SubElement(xmlBox_vectors, 'C')
        xmlCx = ET.SubElement(xmlC, 'x')
        xmlCx.text = str(box_vectors_unitless[2][0])
        xmlCy = ET.SubElement(xmlC, 'y')
        xmlCy.text = str(box_vectors_unitless[2][1])
        xmlCz = ET.SubElement(xmlC, 'z')
        xmlCz.text = str(box_vectors_unitless[2][2])
    else:
        xmlBox_vectors.text = ''
    return

def deserialize_box_vectors(xmlBox_vectors):
    '''
    Takes xml representing box vectors for an openMM or parmed system and
    creates a simtk.unit.Quantity object for it.
    '''
    if xmlBox_vectors.text is not None:
        xmlA = xmlBox_vectors.find('A')
        xmlAx = float(xmlA.find('x').text)
        xmlAy = float(xmlA.find('y').text)
        xmlAz = float(xmlA.find('z').text)
        xmlB = xmlBox_vectors.find('B')
        xmlBx = float(xmlB.find('x').text)
        xmlBy = float(xmlB.find('y').text)
        xmlBz = float(xmlB.find('z').text)
        xmlC = xmlBox_vectors.find('C')
        xmlCx = float(xmlC.find('x').text)
        xmlCy = float(xmlC.find('y').text)
        xmlCz = float(xmlC.find('z').text)
        box_vectors = unit.Quantity([[xmlAx, xmlAy, xmlAz], 
                                     [xmlBx, xmlBy, xmlBz],
                                     [xmlCx, xmlCy, xmlCz]], 
                                     unit=unit.nanometer)
    else:
        box_vectors = None
    return box_vectors

class AmberConfig:
    def __init__(self):
        self.type = "amber"
        self.prmtop_filename = None
        self.inpcrd_filename = None
        self.load_box_vecs_from_coords_file = True
        return
    
    def serialize(self, root):
        assign_tag(root, "type", self.type)
        assign_tag(root, "prmtop_filename", self.prmtop_filename)
        assign_tag(root, "inpcrd_filename", self.inpcrd_filename)
        assign_tag(root, "load_box_vectors_from_coordinates_file", 
                   self.load_box_vecs_from_coords_file)
        return
    
class CharmmConfig:
    def __init__(self):
        self.type = "charmm"
        self.psf_filename = None
        self.pdb_filename = None
        self.params_filenames = []
        return
    
    def serialize(self, root):
        assign_tag(root, "type", self.type)
        assign_tag(root, "psf_filename", self.psf_filename)
        assign_tag(root, "pdb_filename", self.pdb_filename)
        xmlParams = ET.SubElement(root, 'params_filenames')
        for params_filename in self.params_filenames:
            assign_tag(xmlParams, "params_filename", params_filename)
        return

class ForceFieldConfig:
    def __init__(self):
        self.type = "forcefield"
        self.forcefield_filenames = []
        self.pdb_filename = None
        return
    
    def serialize(self, root):
        assign_tag(root, "type", self.type)
        assign_tag(root, "pdb_filename", self.pdb_filename)
        xmlParams = ET.SubElement(root, 'forcefield_filenames')
        for forcefield_filename in self.forcefield_filenames:
            assign_tag(xmlParams, "forcefield_filename", forcefield_filename)
        return
    
class GromacsConfig:
    def __init__(self):
        self.type = "gromacs"
        self.gro_filename = None
        self.top_filename = None
        self.include_dir = None
        return
    
    def serialize(self, root):
        assign_tag(root, "type", self.type)
        assign_tag(root, "gro_filename", self.gro_filename)
        assign_tag(root, "top_filename", self.top_filename)
        assign_tag(root, "include_dir", self.include_dir)
        return

class Config:
    def __init__(self):
        # set all the default values
        
        # system input files: Amber, Charmm, or other
        self.system_files_config = AmberConfig()
        self.box_vectors = None
        self.output_directory = "output/"
        self.overwrite_output = False
        self.chunk_size = 1000
        
        # OpenMM system parameters
        self.nonbonded_method = "PME"
        self.nonbonded_cutoff = 0.9*unit.nanometer
        self.constraints = "HBonds"
        
        # OpenMM integrator parameters
        self.integrator_type = 'langevin'
        self.friction_coefficient = 1.0 / unit.picosecond
        self.target_temperature = 298.15 * unit.kelvin
        self.random_seed = -1
        self.dt = 2.0 * unit.femtosecond
        
        # OpenMM barostat parameters
        self.use_barostat = False
        self.barostat_target_pressure = 1.0 * unit.bar
        self.barostat_target_temperature = 298.15 * unit.kelvin
        self.barostat_frequency = 25
        
        # OpenMM simulation parameters
        self.run_minimization = True
        self.initial_temperature = 298.15 * unit.kelvin

        # OpenMM reporter parameters
        self.energy_reporter_frequency = 5000
        self.coordinates_reporter_frequency = 5000
        self.coordinates_reporter_file_type = "DCD"
        
        # GaMD integrator parameters
        self.gamd_bound = 'lower'
        self.total_simulation_length = 30000
        self.total_boost = True
        self.total_boost_sigma0 = 6.0 * unit.kilocalories_per_mole
        self.dihedral_boost = False
        self.dihedral_boost_sigma0 = 6.0 * unit.kilocalories_per_mole
        self.num_steps_conventional_md = 10000
        self.num_steps_conventional_md_prep = 2000
        self.num_steps_per_averaging = 500
        self.num_steps_gamd_equilibration = 10000
        self.num_steps_gamd_equilibration_prep = 2000
        
        # backup checkpoints
        self.restart_checkpoint_filename = "gamd.backup"
        self.restart_checkpoint_frequency = 1000
        
        self.dihedral_group=1
    
    def serialize(self, filename):
        root = ET.Element('gamd')
        xml_sys_files = ET.SubElement(root, "system_files")
        self.system_files_config.serialize(xml_sys_files)
        serialize_box_vectors(self.box_vectors, root, 
                              tagName='box_vectors')
        assign_tag(root, "output_directory", self.output_directory)
        assign_tag(root, "overwrite_output", self.overwrite_output)
        assign_tag(root, "chunk_size", self.chunk_size)
        assign_tag(root, "nonbonded_method", self.nonbonded_method)
        assign_tag(root, "nonbonded_cutoff", 
                   self.nonbonded_cutoff.value_in_unit(unit.nanometer))
        assign_tag(root, "constraints", self.constraints)
        
        assign_tag(root, "integrator_type", self.integrator_type)
        assign_tag(root, "friction_coefficient", 
                   self.friction_coefficient.value_in_unit(unit.picosecond**-1))
        assign_tag(root, "target_temperature", 
                   self.target_temperature.value_in_unit(unit.kelvin))
        assign_tag(root, "random_seed", self.random_seed)
        assign_tag(root, "dt", 
                   self.dt.value_in_unit(unit.picosecond))
        
        assign_tag(root, "use_barostat", self.use_barostat)
        assign_tag(root, "barostat_target_pressure", 
                   self.barostat_target_pressure.value_in_unit(unit.bar))
        assign_tag(root, "barostat_target_temperature", 
                   self.barostat_target_temperature.value_in_unit(unit.kelvin))
        assign_tag(root, "barostat_frequency", self.barostat_frequency)
        
        assign_tag(root, "run_minimization", self.run_minimization)
        assign_tag(root, "initial_temperature", 
                   self.initial_temperature.value_in_unit(unit.kelvin))
        
        assign_tag(root, "energy_reporter_frequency", 
                   self.energy_reporter_frequency)
        assign_tag(root, "coordinates_reporter_frequency", 
                   self.coordinates_reporter_frequency)
        assign_tag(root, "coordinates_reporter_file_type", 
                   self.coordinates_reporter_file_type)
        
        assign_tag(root, "gamd_bound", self.gamd_bound)
        assign_tag(root, "total_simulation_length", 
                   self.total_simulation_length)
        assign_tag(root, "total_boost", self.total_boost)
        assign_tag(root, "total_boost_sigma0", 
                   self.total_boost_sigma0.value_in_unit(
                       unit.kilocalories_per_mole))
        assign_tag(root, "dihedral_boost", self.dihedral_boost)
        assign_tag(root, "dihedral_boost_sigma0", 
                   self.dihedral_boost_sigma0.value_in_unit(
                       unit.kilocalories_per_mole))
        assign_tag(root, "num_steps_conventional_md", 
                   self.num_steps_conventional_md)
        assign_tag(root, "num_steps_conventional_md_prep", 
                   self.num_steps_conventional_md_prep)
        assign_tag(root, "num_steps_per_averaging", 
                   self.num_steps_per_averaging)
        assign_tag(root, "num_steps_gamd_equilibration", 
                   self.num_steps_gamd_equilibration)
        assign_tag(root, "num_steps_gamd_equilibration_prep", 
                   self.num_steps_gamd_equilibration_prep)
        
        assign_tag(root, "restart_checkpoint_filename", 
                   self.restart_checkpoint_filename)
        assign_tag(root, "restart_checkpoint_frequency", 
                   self.restart_checkpoint_frequency)
        
        assign_tag(root, "dihedral_group", self.dihedral_group)
        
        xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(
            indent="   ")
        our_file=open(filename, 'w')
        our_file.write(xmlstr)
        our_file.close()
        

if __name__ == "__main__":
    myconfig = Config()
    myconfig.serialize("/tmp/gamdconfig.xml")