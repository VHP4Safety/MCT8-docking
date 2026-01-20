"""
Tests for Flask application endpoints.
"""

import pytest
from app import allowed_file, validate_docking_params


class TestAllowedFile:
    """Tests for allowed_file function."""

    def test_allowed_extensions(self):
        """Test that allowed extensions pass."""
        assert allowed_file("molecule.sdf") is True
        assert allowed_file("smiles.smi") is True
        assert allowed_file("data.csv") is True
        assert allowed_file("input.txt") is True

    def test_disallowed_extensions(self):
        """Test that disallowed extensions fail."""
        assert allowed_file("script.py") is False
        assert allowed_file("data.json") is False
        assert allowed_file("image.png") is False
        assert allowed_file("document.pdf") is False

    def test_no_extension(self):
        """Test files without extension."""
        assert allowed_file("noextension") is False

    def test_case_insensitive(self):
        """Test that extension check is case-insensitive."""
        assert allowed_file("MOLECULE.SDF") is True
        assert allowed_file("Data.CSV") is True
        assert allowed_file("file.TXT") is True

    def test_double_extension(self):
        """Test files with double extensions."""
        assert allowed_file("file.tar.sdf") is True
        assert allowed_file("file.sdf.py") is False


class TestValidateDockingParams:
    """Tests for validate_docking_params function."""

    def test_valid_params(self):
        """Test that valid parameters pass validation."""
        params = {
            'num_modes': 3,
            'exhaustiveness': 8,
            'autobox_add': 6.0,
            'cnn': 'none'
        }
        errors = validate_docking_params(params)
        assert len(errors) == 0

    def test_num_modes_bounds(self):
        """Test num_modes parameter bounds."""
        # Too low
        errors = validate_docking_params({'num_modes': 0})
        assert len(errors) == 1
        assert "num_modes" in errors[0]

        # Too high
        errors = validate_docking_params({'num_modes': 51})
        assert len(errors) == 1
        assert "num_modes" in errors[0]

        # Valid boundaries
        errors = validate_docking_params({'num_modes': 1})
        assert len(errors) == 0
        errors = validate_docking_params({'num_modes': 50})
        assert len(errors) == 0

    def test_exhaustiveness_bounds(self):
        """Test exhaustiveness parameter bounds."""
        # Too low
        errors = validate_docking_params({'exhaustiveness': 1})
        assert len(errors) == 1
        assert "exhaustiveness" in errors[0]

        # Too high
        errors = validate_docking_params({'exhaustiveness': 129})
        assert len(errors) == 1
        assert "exhaustiveness" in errors[0]

        # Valid boundaries
        errors = validate_docking_params({'exhaustiveness': 2})
        assert len(errors) == 0
        errors = validate_docking_params({'exhaustiveness': 128})
        assert len(errors) == 0

    def test_autobox_add_bounds(self):
        """Test autobox_add parameter bounds."""
        # Too low
        errors = validate_docking_params({'autobox_add': -1})
        assert len(errors) == 1
        assert "autobox_add" in errors[0]

        # Too high
        errors = validate_docking_params({'autobox_add': 11})
        assert len(errors) == 1
        assert "autobox_add" in errors[0]

        # Valid boundaries
        errors = validate_docking_params({'autobox_add': 0})
        assert len(errors) == 0
        errors = validate_docking_params({'autobox_add': 10})
        assert len(errors) == 0

    def test_cnn_values(self):
        """Test cnn parameter allowed values."""
        # Valid values
        for cnn in ['none', 'fast', 'default']:
            errors = validate_docking_params({'cnn': cnn})
            assert len(errors) == 0

        # Invalid value
        errors = validate_docking_params({'cnn': 'invalid'})
        assert len(errors) == 1
        assert "cnn" in errors[0]

    def test_multiple_errors(self):
        """Test that multiple invalid params produce multiple errors."""
        params = {
            'num_modes': 0,
            'exhaustiveness': 1,
            'autobox_add': -5,
            'cnn': 'invalid'
        }
        errors = validate_docking_params(params)
        assert len(errors) == 4

    def test_empty_params(self):
        """Test empty parameters dictionary."""
        errors = validate_docking_params({})
        assert len(errors) == 0


class TestHomeRoute:
    """Tests for home route."""

    def test_home_returns_200(self, client):
        """Test that home page loads successfully."""
        response = client.get('/')
        assert response.status_code == 200

    def test_home_contains_form(self, client):
        """Test that home page contains the docking form."""
        response = client.get('/')
        assert b'smiles' in response.data or b'SMILES' in response.data


class TestApiEndpoints:
    """Tests for API endpoints."""

    def test_api_missing_smiles(self, client):
        """Test API returns 400 when SMILES is missing."""
        response = client.post('/api',
                              json={},
                              content_type='application/json')
        assert response.status_code == 400
        data = response.get_json()
        assert 'error' in data

    def test_api_smiles_not_list(self, client):
        """Test API returns 400 when SMILES is not a list."""
        response = client.post('/api',
                              json={'smiles': 'CCO'},
                              content_type='application/json')
        assert response.status_code == 400
        data = response.get_json()
        assert 'error' in data

    def test_api_invalid_smiles(self, client):
        """Test API returns 400 when all SMILES are invalid."""
        response = client.post('/api',
                              json={'smiles': ['INVALID', 'NOT_SMILES']},
                              content_type='application/json')
        assert response.status_code == 400
        data = response.get_json()
        assert 'error' in data
        assert 'invalid_smiles' in data

    def test_api_invalid_params(self, client):
        """Test API returns 400 for invalid parameters."""
        response = client.post('/api',
                              json={
                                  'smiles': ['CCO'],
                                  'params': {'num_modes': 100}
                              },
                              content_type='application/json')
        assert response.status_code == 400
        data = response.get_json()
        assert 'error' in data
        assert 'num_modes' in data['error']

    def test_api_too_many_smiles(self, client):
        """Test API returns 400 when too many SMILES provided."""
        smiles_list = ['CCO'] * 101  # MAX_SMILES_INPUT is 100
        response = client.post('/api',
                              json={'smiles': smiles_list},
                              content_type='application/json')
        assert response.status_code == 400
        data = response.get_json()
        assert 'error' in data
        assert '100' in data['error']


class TestReceptorEndpoint:
    """Tests for receptor API endpoint."""

    def test_receptor_endpoint_without_file(self, client):
        """Test receptor endpoint when PDB file doesn't exist."""
        # This test may pass or fail depending on whether
        # the receptor file exists in the test environment
        response = client.get('/api/receptor')
        # Either 200 (file exists) or 404 (file missing)
        assert response.status_code in [200, 404]


class TestBindingSiteEndpoint:
    """Tests for binding site API endpoint."""

    def test_binding_site_endpoint_without_file(self, client):
        """Test binding site endpoint when PDB file doesn't exist."""
        response = client.get('/api/binding_site')
        # Either 200 (file exists) or 404 (file missing)
        assert response.status_code in [200, 404]


class TestDownloadEndpoint:
    """Tests for download endpoint."""

    def test_download_without_results(self, client):
        """Test download endpoint when no results exist."""
        response = client.get('/download/sdf')
        # Should redirect to home if no results
        assert response.status_code == 302
