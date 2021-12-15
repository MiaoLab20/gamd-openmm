"""
test_config.py

Test the config.py module.
"""

import os

import pytest

from gamd import config

TEST_DIRECTORY = os.path.dirname(__file__)

def test_default_config_write(tmp_path, default_config):
    """
    Test the ability to create and write a Config object.
    """
    default_config_path = os.path.join(tmp_path, "default_config.xml")
    default_config.serialize(default_config_path)
    assert os.path.exists(default_config_path)
    return