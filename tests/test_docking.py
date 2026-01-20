"""
Unit tests for docking.py module.
"""

import pytest
import docking


class TestValidateSmiles:
    """Tests for validate_smiles function."""

    def test_valid_smiles(self, valid_smiles):
        """Test that valid SMILES are recognized."""
        result = docking.validate_smiles(valid_smiles)
        assert len(result['valid']) == 3
        assert len(result['invalid']) == 0
        assert "CCO" in result['valid']
        assert "c1ccccc1" in result['valid']

    def test_invalid_smiles(self, invalid_smiles):
        """Test that invalid SMILES are rejected."""
        result = docking.validate_smiles(invalid_smiles)
        assert len(result['valid']) == 0
        assert len(result['invalid']) == 3
        assert "INVALID" in result['invalid']

    def test_mixed_smiles(self, valid_smiles, invalid_smiles):
        """Test mixed valid and invalid SMILES."""
        mixed = valid_smiles + invalid_smiles
        result = docking.validate_smiles(mixed)
        assert len(result['valid']) == 3
        assert len(result['invalid']) == 3

    def test_empty_list(self):
        """Test empty input list."""
        result = docking.validate_smiles([])
        assert len(result['valid']) == 0
        assert len(result['invalid']) == 0

    def test_whitespace_handling(self):
        """Test that whitespace is stripped from SMILES."""
        result = docking.validate_smiles(["  CCO  ", "\tc1ccccc1\n"])
        assert len(result['valid']) == 2
        assert "CCO" in result['valid']
        assert "c1ccccc1" in result['valid']


class TestAssessInhibition:
    """Tests for assess_inhibition function."""

    def test_likely_inhibitor(self, sample_affinities):
        """Test classification of likely inhibitors."""
        result = docking.assess_inhibition(sample_affinities['likely_inhibitor'])
        assert result['category'] == 'Likely Inhibitor'
        assert result['color'] == '#E6007E'

    def test_possible_inhibitor(self, sample_affinities):
        """Test classification of possible inhibitors."""
        result = docking.assess_inhibition(sample_affinities['possible_inhibitor'])
        assert result['category'] == 'Possible Inhibitor'
        assert result['color'] == '#FF9500'

    def test_unlikely_inhibitor(self, sample_affinities):
        """Test classification of unlikely inhibitors."""
        result = docking.assess_inhibition(sample_affinities['unlikely_inhibitor'])
        assert result['category'] == 'Unlikely Inhibitor'
        assert result['color'] == '#4CAF50'

    def test_boundary_at_minus_nine(self, sample_affinities):
        """Test boundary at -9.0 kcal/mol (should be Possible)."""
        result = docking.assess_inhibition(sample_affinities['boundary_likely'])
        # -9.0 is NOT < -9.0, so it's Possible Inhibitor
        assert result['category'] == 'Possible Inhibitor'

    def test_boundary_at_minus_eight(self, sample_affinities):
        """Test boundary at -8.0 kcal/mol (should be Unlikely)."""
        result = docking.assess_inhibition(sample_affinities['boundary_possible'])
        # -8.0 is NOT < -8.0, so it's Unlikely Inhibitor
        assert result['category'] == 'Unlikely Inhibitor'

    def test_none_affinity(self):
        """Test handling of None affinity."""
        result = docking.assess_inhibition(None)
        assert result['category'] == 'Unknown'
        assert result['color'] == '#999999'

    def test_zero_affinity(self):
        """Test zero affinity (should be unlikely)."""
        result = docking.assess_inhibition(0)
        assert result['category'] == 'Unlikely Inhibitor'


class TestSmilesToImage:
    """Tests for smiles_to_image function."""

    def test_valid_smiles_returns_data_uri(self):
        """Test that valid SMILES returns a data URI."""
        result = docking.smiles_to_image("CCO")
        assert result is not None
        assert result.startswith("data:image/png;base64,")

    def test_invalid_smiles_returns_none(self):
        """Test that invalid SMILES returns None."""
        result = docking.smiles_to_image("INVALID")
        assert result is None

    def test_custom_image_size(self):
        """Test image generation with custom size."""
        result = docking.smiles_to_image("CCO", img_size=(100, 100))
        assert result is not None
        assert result.startswith("data:image/png;base64,")
