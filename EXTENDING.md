# Extending the MCT8 Docking Application

This guide explains how to extend the application to support new receptors, customize scoring thresholds, and add new features.

## Table of Contents

1. [Adding a New Receptor](#adding-a-new-receptor)
2. [Configuring Binding Sites](#configuring-binding-sites)
3. [Customizing Scoring and Assessment](#customizing-scoring-and-assessment)
4. [Testing Your Changes](#testing-your-changes)

---

## Adding a New Receptor

### Step 1: Prepare the PDB File

Obtain your receptor structure from one of these sources:
- **RCSB PDB**: https://www.rcsb.org/ (experimental structures)
- **AlphaFold DB**: https://alphafold.ebi.ac.uk/ (predicted structures)
- **Homology modeling**: Tools like SWISS-MODEL or Modeller

Ensure the PDB file:
- Contains only protein atoms (remove water, ions, co-crystallized ligands)
- Has correct protonation states for the pH of interest
- Is energy-minimized if from homology modeling

### Step 2: Add Receptor to docking.py

Open `docking.py` and add your receptor as a constant. The MCT8 receptor is defined starting at line 74:

```python
# Add your receptor constant (similar to MCT8_RECEPTOR_PDB)
MY_RECEPTOR_PDB = """ATOM      1  N   MET A   1      ...
ATOM      2  CA  MET A   1      ...
...
END
"""
```

### Step 3: Update create_pdb_files()

Modify the `create_pdb_files()` function to write your new receptor:

```python
def create_pdb_files():
    """Create PDB files for docking."""
    os.makedirs('data', exist_ok=True)

    # Write MCT8 receptor (existing)
    with open('data/mct8_receptor_full.pdb', 'w') as f:
        f.write(MCT8_RECEPTOR_PDB)

    # Write your new receptor
    with open('data/my_receptor.pdb', 'w') as f:
        f.write(MY_RECEPTOR_PDB)

    # ... rest of function
```

### Step 4: Add Receptor Selection to UI (Optional)

If you want users to choose between receptors, update `templates/index.html` to add a receptor dropdown and modify `app.py` to handle the selection.

---

## Configuring Binding Sites

The binding site defines where Gnina searches for ligand binding poses. It's specified as a set of coordinates in PDB HETATM format.

### Understanding Binding Site Format

The binding site in `docking.py` (lines 33-55) uses HETATM records:

```
HETATM    1  C   SIT E   1      -0.639  23.495  -2.246  1.00  0.00           C
```

Fields:
- `HETATM`: Record type for heteroatoms
- `1`: Atom serial number
- `C`: Atom name
- `SIT`: Residue name (arbitrary, "SIT" for site)
- `E`: Chain identifier
- `1`: Residue number
- `-0.639  23.495  -2.246`: X, Y, Z coordinates (Ångstroms)
- `1.00  0.00`: Occupancy and temperature factor
- `C`: Element symbol

### Extracting Binding Site from Known Ligand

If you have a crystal structure with a bound ligand:

1. **Using PyMOL**:
   ```python
   # In PyMOL
   load receptor_with_ligand.pdb
   select ligand, resn LIG  # Replace LIG with ligand residue name
   save binding_site.pdb, ligand
   ```

2. **Using RDKit/Python**:
   ```python
   from rdkit import Chem

   # Load structure and extract ligand coordinates
   mol = Chem.MolFromPDBFile('receptor_with_ligand.pdb')
   # Get coordinates of ligand atoms
   ```

3. **From literature**: Use published binding site residue information to define the search box center.

### Defining Binding Site Manually

If no bound ligand is available, define the binding site from:
- Known active site residues from literature
- Molecular dynamics simulations
- Cavity detection tools (fpocket, SiteMap)

Create coordinates representing the approximate center and extent of the binding pocket.

---

## Customizing Scoring and Assessment

### Threshold Constants

Risk assessment thresholds are defined in `docking.py` at lines 29-30:

```python
STRONG_INHIBITOR_THRESHOLD = -9.0  # kcal/mol
WEAK_INHIBITOR_THRESHOLD = -8.0    # kcal/mol
```

Adjust these based on:
- Experimental validation data for your receptor
- Literature benchmarks
- Desired sensitivity/specificity tradeoff

### Modifying assess_inhibition()

The `assess_inhibition()` function (lines 382-419 in `docking.py`) categorizes compounds:

```python
def assess_inhibition(affinity):
    """
    Assess inhibition risk based on binding affinity.

    Args:
        affinity: Binding affinity in kcal/mol

    Returns:
        dict with category, color, and description
    """
    if affinity < STRONG_INHIBITOR_THRESHOLD:
        return {
            'category': 'Likely Inhibitor',
            'color': '#E6007E',  # VHP4Safety magenta
            'description': f'High inhibition risk (< {STRONG_INHIBITOR_THRESHOLD} kcal/mol)'
        }
    elif affinity < WEAK_INHIBITOR_THRESHOLD:
        return {
            'category': 'Possible Inhibitor',
            'color': '#FF9500',  # Orange
            'description': f'Moderate inhibition risk ({WEAK_INHIBITOR_THRESHOLD} to {STRONG_INHIBITOR_THRESHOLD} kcal/mol)'
        }
    else:
        return {
            'category': 'Unlikely Inhibitor',
            'color': '#4CAF50',  # Green
            'description': f'Low inhibition risk (> {WEAK_INHIBITOR_THRESHOLD} kcal/mol)'
        }
```

### Updating Colors in CSS

Risk assessment colors are also defined in `static/styles.css` (around line 380):

```css
.badge-likely {
    background-color: #E6007E;  /* Magenta - high risk */
}

.badge-possible {
    background-color: #FF9500;  /* Orange - moderate risk */
}

.badge-unlikely {
    background-color: #4CAF50;  /* Green - low risk */
}
```

Ensure CSS colors match the Python assessment colors for consistency.

---

## Testing Your Changes

### Running Unit Tests

The test suite in `tests/` covers core functionality without requiring Gnina:

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=. --cov-report=term-missing

# Run specific test file
pytest tests/test_docking.py -v
```

### Adding New Tests

Add tests for new functionality in `tests/test_docking.py`:

```python
def test_my_new_receptor_validation():
    """Test that new receptor PDB is valid."""
    from docking import MY_RECEPTOR_PDB

    # Check PDB has content
    assert len(MY_RECEPTOR_PDB) > 0

    # Check has ATOM records
    assert 'ATOM' in MY_RECEPTOR_PDB

    # Check ends with END
    assert MY_RECEPTOR_PDB.strip().endswith('END')


def test_custom_threshold():
    """Test custom inhibition thresholds."""
    from docking import assess_inhibition

    # Test boundary conditions
    result = assess_inhibition(-9.1)
    assert result['category'] == 'Likely Inhibitor'

    result = assess_inhibition(-8.5)
    assert result['category'] == 'Possible Inhibitor'

    result = assess_inhibition(-7.0)
    assert result['category'] == 'Unlikely Inhibitor'
```

### Manual Testing

1. **Start the Flask application**:
   ```bash
   python app.py
   ```

2. **Test with sample SMILES**:
   - Simple: `CCO` (ethanol)
   - Complex: `c1ccccc1` (benzene)
   - Default: Silychristin (pre-loaded)

3. **Verify 3D viewer**:
   - Click a result row to open viewer
   - Check receptor loads correctly
   - Verify binding site displays
   - Test multi-pose comparison

4. **Check browser console**:
   - Open DevTools (F12)
   - Look for JavaScript errors
   - Verify API responses

### Docker Testing

After changes, rebuild and test the Docker image:

```bash
# Build image
docker build -t mct8-docking:test .

# Run container
docker run -d -p 5000:5000 --name test-container mct8-docking:test

# Test endpoints
curl http://localhost:5000/api/receptor | head -c 200
curl http://localhost:5000/api/binding_site

# Cleanup
docker rm -f test-container
```

---

## Key Extension Points Reference

| File | Lines | Purpose |
|------|-------|---------|
| `docking.py` | 29-30 | Threshold constants |
| `docking.py` | 33-55 | BINDING_SITE_PDB definition |
| `docking.py` | 74-10500 | MCT8_RECEPTOR_PDB definition |
| `docking.py` | 382-419 | `assess_inhibition()` function |
| `docking.py` | 280-330 | `run_docking()` Gnina parameters |
| `app.py` | 38-70 | Parameter validation |
| `app.py` | 150-200 | `/dock` route handler |
| `app.py` | 250-300 | `/api` route handler |
| `templates/index.html` | 100-200 | Input form controls |
| `templates/index.html` | 400-500 | 3D viewer JavaScript |
| `static/styles.css` | 370-400 | Risk assessment badge colors |

---

## Common Extension Scenarios

### Scenario 1: Add Multiple Receptor Support

1. Add receptor constants to `docking.py`
2. Create receptor dropdown in `templates/index.html`
3. Modify `/dock` and `/api` routes to accept receptor parameter
4. Update `run_docking()` to use selected receptor file

### Scenario 2: Add Custom Scoring Function

1. Implement scoring function in `docking.py`
2. Call from `process_results()` after parsing SDF
3. Add score to result dictionary
4. Display in results table and 3D viewer info panel

### Scenario 3: Add New Output Format

1. Add format option to API parameter validation in `app.py`
2. Implement format conversion in `/api` route
3. Update API documentation in README

---

## Getting Help

- **GitHub Issues**: Report bugs or request features
- **Gnina Documentation**: https://github.com/gnina/gnina
- **RDKit Documentation**: https://www.rdkit.org/docs/
- **VHP4Safety**: https://www.vhp4safety.nl/
