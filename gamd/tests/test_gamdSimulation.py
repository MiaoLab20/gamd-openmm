"""
test_gamdSimulation.py

Test the gamdSimulation.py module.
"""

import os

import pytest
import openmm.unit as unit

from gamd import parser
from gamd import gamdSimulation

TEST_DIRECTORY = os.path.dirname(__file__)
ROOT_DIRECTORY = os.path.join(TEST_DIRECTORY, "..")

def test_amber_alanine_dipeptide_object(tmp_path):
    """
    
    """
    os.chdir(ROOT_DIRECTORY)
    input_file = os.path.join(TEST_DIRECTORY, "data/dip_amber.xml")
    parserFactory = parser.ParserFactory()
    config = parserFactory.parse_file(input_file, "xml")
    gamdSimulationFactory = gamdSimulation.GamdSimulationFactory()
    gamdSim = gamdSimulationFactory.createGamdSimulation(
        config, "reference", "")