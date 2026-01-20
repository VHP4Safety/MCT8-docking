"""
Pytest fixtures for MCT8 docking tests.
"""

import os
import sys
import pytest

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app import app as flask_app


@pytest.fixture
def app():
    """Create application for testing."""
    flask_app.config.update({
        "TESTING": True,
    })
    yield flask_app


@pytest.fixture
def client(app):
    """Create test client."""
    return app.test_client()


@pytest.fixture
def valid_smiles():
    """Return list of valid SMILES strings."""
    return [
        "CCO",        # Ethanol
        "c1ccccc1",   # Benzene
        "CC(=O)O",    # Acetic acid
    ]


@pytest.fixture
def invalid_smiles():
    """Return list of invalid SMILES strings."""
    return [
        "INVALID",
        "NOT_A_SMILES",
        "XYZ123",
    ]


@pytest.fixture
def sample_affinities():
    """Return sample affinity values for testing thresholds."""
    return {
        "likely_inhibitor": -10.5,    # < -9.0
        "possible_inhibitor": -8.5,    # -8.0 to -9.0
        "unlikely_inhibitor": -6.0,    # > -8.0
        "boundary_likely": -9.0,       # Exactly at threshold
        "boundary_possible": -8.0,     # Exactly at threshold
    }
