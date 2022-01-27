"""
test_parser.py

Test the parser.py module.
"""

import os

import pytest
import openmm.unit as unit

from gamd import parser

TEST_DIRECTORY = os.path.dirname(__file__)

def test_amber_alanine_dipeptide_read_rewrite(tmp_path):
    """
    
    """
    input_file = os.path.join(TEST_DIRECTORY, "data/dip_amber.xml")
    output_file = os.path.join(tmp_path, "dip_amber_rewrite.xml")
    myparser = parser.XmlParser()
    myparser.parse_file(input_file)
    myparser.config.serialize(output_file)
    assert os.path.exists(output_file)
    myparser2 = parser.XmlParser()
    myparser2.parse_file(output_file)
    config = myparser2.config
    check_config_from_test_input_xml(config)
    check_config_from_test_amber_input_xml(config)
    

def check_config_from_test_amber_input_xml(config):
    """
    
    """
    assert config.input_files.amber.topology == "data/dip.top"
    assert config.input_files.amber.coordinates == "data/md-4ns.rst7"

def check_config_from_test_input_xml(config):
    """
    
    """
    assert config.temperature == 298.15 * unit.kelvin
    assert config.system.nonbonded_method == "pme"
    assert config.system.nonbonded_cutoff == 1.0 * unit.nanometers
    assert config.system.constraints == "hbonds"
    assert config.barostat.pressure == 1.0 * unit.bar
    assert config.barostat.frequency == 25
    assert config.run_minimization == True
    assert config.integrator.algorithm == "langevin"
    assert config.integrator.boost_type == "lower-dual"
    assert config.integrator.sigma0.primary == 6.0 * unit.kilocalories_per_mole
    assert config.integrator.sigma0.secondary == 6.0 * unit.kilocalories_per_mole
    assert config.integrator.random_seed == 0
    assert config.integrator.dt == 0.002 * unit.picoseconds
    assert config.integrator.friction_coefficient == 1.0 * unit.picoseconds ** -1
    assert config.integrator.number_of_steps.conventional_md_prep == 2000
    assert config.integrator.number_of_steps.conventional_md == 10000
    assert config.integrator.number_of_steps.gamd_equilibration_prep == 2000
    assert config.integrator.number_of_steps.gamd_equilibration == 20000
    assert config.integrator.number_of_steps.gamd_production == 30000
    assert config.integrator.number_of_steps.averaging_window_interval == 50
    assert config.outputs.directory == "output/"
    assert config.outputs.overwrite_output == True
    assert config.outputs.reporting.energy_interval == 500
    assert config.outputs.reporting.coordinates_file_type == "dcd"
    assert config.outputs.reporting.coordinates_interval == 500
    assert config.outputs.reporting.restart_checkpoint_interval == 50000
    assert config.outputs.reporting.statistics_interval == 500
    return

def test_amber_alanine_dipeptide_read_check(tmp_path):
    """
    
    """
    input_file = os.path.join(TEST_DIRECTORY, "data/dip_amber.xml")
    myparser = parser.XmlParser()
    myparser.parse_file(input_file)
    config = myparser.config
    check_config_from_test_input_xml(config)
    check_config_from_test_amber_input_xml(config)
    
    