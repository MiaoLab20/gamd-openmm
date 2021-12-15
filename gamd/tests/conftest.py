"""
conftest.py

configurations for all tests
"""

import os

import pytest

from gamd import config

TEST_DIRECTORY = os.path.dirname(__file__)

@pytest.fixture(scope="session")
def default_config():
    """
    Create a Config object with all the default values.
    """
    myconfig = config.Config()
    return myconfig

@pytest.fixture(scope="session")
def amber_alanine_dipeptide_config():
    """
    Create a Config object for alanine dipeptide using the AMBER forcefield.
    """
    
    return 