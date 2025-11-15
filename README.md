# MCT8 Docking - VHP4Safety

Molecular docking web application for predicting inhibitor binding to Monocarboxylate Transporter 8 (MCT8), part of the Virtual Human Platform for Safety (VHP4Safety) thyroid case study.

## Overview

MCT8 (Monocarboxylate Transporter 8) is critical for thyroid hormone transport during early embryonic development. Inhibition of MCT8 can lead to developmental disorders and other adverse outcomes. This tool helps assess the inhibition risk of chemical compounds using neural network-enhanced molecular docking with Gnina.

## Features

- **Web Interface**: User-friendly form for submitting ligands and configuring docking parameters
- **Multiple Input Methods**: SMILES strings, SDF files, or CSV uploads
- **Automated Ligand Preparation**: Hydrogen addition, 3D structure generation, and geometry optimization
- **Neural Network Scoring**: Gnina's CNN-enhanced scoring for improved accuracy
- **Interactive Results**: DataTables with sorting, filtering, and export (CSV, Excel, PDF)
- **Risk Assessment**: Automatic categorization of compounds as strong/moderate/weak inhibitors
- **RESTful API**: Programmatic access for batch processing
- **PDF Reports**: Generate publication-ready docking reports

## Quick Start

### Docker Deployment (Recommended)

```bash
# Build the Docker image
docker build -t mct8-docking .

# Run the container
docker run -d -p 5000:5000 --name mct8-docking mct8-docking

# Access the application
open http://localhost:5000
```

### Local Development

```bash
# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Download Gnina (done automatically on first run)
python -c "import docking; docking.setup_environment()"

# Run the Flask application
python app.py

# Access at http://localhost:5000
```

## Usage

### Web Interface

1. **Input Ligands**:
   - Enter SMILES strings (comma-separated)
   - Or upload SDF, SMI, or CSV file (with 'SMILES' column)
   - Default: Silychristin (known MCT8 inhibitor)

2. **Configure Parameters**:
   - **Number of Poses**: 1-50 (default: 3)
   - **Exhaustiveness**: 2-128 (default: 8) - higher = more thorough search
   - **Search Box**: 0-10 Å (default: 6) - expansion around binding site
   - **CNN Scoring**: fast/default/none

3. **Run Docking**:
   - Click "Run Docking" to execute simulation
   - Or "Generate PDF Report" for downloadable report

4. **Interpret Results**:
   - **Strong Inhibitors**: Affinity < -9.0 kcal/mol (High risk)
   - **Moderate Inhibitors**: -8.0 to -9.0 kcal/mol (Moderate risk)
   - **Weak/Non-Inhibitors**: > -8.0 kcal/mol (Low risk)

### REST API

#### Endpoint: `POST /api`

**Request:**
```json
{
  "smiles": [
    "COC1=C(C=CC(=C1)[C@H]2[C@@H]...",
    "c1ccccc1"
  ],
  "params": {
    "num_modes": 3,
    "exhaustiveness": 8,
    "autobox_add": 6.0,
    "cnn": "fast"
  },
  "format": "json"
}
```

**Response:**
```json
{
  "success": true,
  "num_poses": 6,
  "invalid_smiles": [],
  "results": [
    {
      "pose_index": 0,
      "smiles": "COC1=C...",
      "inchikey": "ZQSIJRDFPHDXBS-UHFFFAOYSA-N",
      "affinity": -9.5,
      "cnn_score": 0.842,
      "boltzmann_weight": 0.453,
      "assessment": "Strong Inhibitor"
    }
  ]
}
```

**cURL Example:**
```bash
curl -X POST http://localhost:5000/api \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": ["CCO", "c1ccccc1"],
    "params": {"num_modes": 5},
    "format": "json"
  }'
```

**Output Formats:**
- `"format": "json"` - JSON response (default)
- `"format": "csv"` - CSV file
- `"format": "sdf"` - SDF file with 3D structures

## Project Structure

```
MCT8-docking/
├── app.py                      # Flask application (routes, API)
├── docking.py                  # Core docking logic
├── requirements.txt            # Python dependencies
├── Dockerfile                  # Docker configuration
├── entrypoint.sh              # Container startup script
├── README.md                   # This file
├── .gitignore                 # Git exclusions
├── static/                     # Static assets
│   ├── styles.css             # VHP4Safety styling
│   └── img/
│       ├── logo.png           # VHP4Safety logo
│       └── github.svg         # GitHub icon
├── templates/                  # HTML templates
│   └── index.html             # Main UI
├── data/                       # Static data files
│   ├── mct8_receptor_full.pdb # MCT8 protein structure
│   └── binding_site.pdb       # Binding site (auto-generated)
├── binaries/                   # External binaries
│   └── gnina                   # Gnina executable (auto-downloaded)
└── samples/                    # Sample data
    └── silychristin.smi       # Default test molecule
```

## Technical Details

### Docking Workflow

1. **Ligand Validation**: SMILES strings validated with RDKit
2. **3D Structure Generation**: ETKDGv3 embedding algorithm
3. **Geometry Optimization**: MMFF force field minimization
4. **Molecular Docking**: Gnina with CNN scoring
5. **Result Processing**: Boltzmann weighting, affinity assessment
6. **Visualization**: 2D structure images, 3D pose viewing (future)

### Scoring

- **Affinity**: Binding free energy (kcal/mol) - primary metric
- **CNN Score**: Neural network confidence (0-1)
- **Boltzmann Weight**: Thermodynamic probability of pose

### MCT8 Receptor Model

The MCT8 protein structure is based on computational modeling for the thyroid case study. The binding site coordinates are defined in `data/binding_site.pdb`.

## Dependencies

### Core Libraries
- **Flask**: Web framework
- **RDKit**: Cheminformatics and molecular processing
- **Gnina**: Neural network-enhanced docking
- **pandas**: Data manipulation
- **reportlab**: PDF generation

### System Requirements
- Python 3.10+
- 2GB RAM minimum (4GB recommended for larger batches)
- Docker (for containerized deployment)

## Development

### Running Tests
```bash
# Test docking module
python docking.py

# Test with sample data
python -c "import docking; print(docking.validate_smiles(['CCO', 'INVALID']))"
```

### Adding New Features
1. Create feature branch: `git checkout -b feature/description`
2. Make changes and test thoroughly
3. Commit with descriptive message
4. Push and create pull request

### Code Style
- PEP 8 compliance
- Type hints where applicable
- Docstrings for all functions
- Logging for debugging

## Troubleshooting

### Gnina Download Fails
```bash
# Manually download Gnina
wget https://github.com/gnina/gnina/releases/download/v1.3/gnina
chmod +x gnina
mv gnina binaries/
```

### Docker Build Issues
```bash
# Clear Docker cache
docker system prune -a

# Rebuild without cache
docker build --no-cache -t mct8-docking .
```

### RDKit Import Error
```bash
# Install RDKit conda package (recommended)
conda install -c conda-forge rdkit

# Or use pip (may require system dependencies)
pip install rdkit
```

## Citation

If you use this tool in your research, please cite:

```
VHP4Safety MCT8 Docking Tool
Virtual Human Platform for Safety - Thyroid Case Study
https://github.com/VHP4Safety/mct8-docking
```

## License

This project is part of the VHP4Safety initiative. See LICENSE file for details.

## Acknowledgments

- **Gnina**: Koes et al. - Neural network-enhanced molecular docking
- **RDKit**: Open-source cheminformatics toolkit
- **VHP4Safety**: Virtual Human Platform for Safety project
- **Original Notebook**: MCT8 docking Colab notebook contributors

## Contact

For issues, questions, or contributions:
- GitHub Issues: https://github.com/VHP4Safety/mct8-docking/issues
- VHP4Safety: https://www.vhp4safety.nl/

## Related Links

- [VHP4Safety Project](https://www.vhp4safety.nl/)
- [Gnina Documentation](https://github.com/gnina/gnina)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [MCT8 (SLC16A2) on Wikipedia](https://en.wikipedia.org/wiki/SLC16A2)
