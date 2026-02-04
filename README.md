# MCT8 Docking - VHP4Safety

[![CI](https://github.com/VHP4Safety/MCT8-docking/actions/workflows/ci.yml/badge.svg)](https://github.com/VHP4Safety/MCT8-docking/actions/workflows/ci.yml)
[![Docker Build](https://github.com/VHP4Safety/MCT8-docking/actions/workflows/docker.yml/badge.svg)](https://github.com/VHP4Safety/MCT8-docking/actions/workflows/docker.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://www.docker.com/)

Molecular docking web application for predicting binding affinity to Monocarboxylate Transporter 8 (MCT8), part of the Virtual Human Platform for Safety (VHP4Safety) thyroid case study.

## Overview

MCT8 (Monocarboxylate Transporter 8, also known as SLC16A2) is essential for thyroid hormone transport across cell membranes, particularly critical during early embryonic brain development. MCT8 facilitates the uptake of thyroid hormones (T3 and T4) into developing neurons, which is vital for proper neurological development in the first trimester.

Inhibition of MCT8 can disrupt thyroid hormone availability in the developing brain, potentially leading to developmental disorders and adverse neurodevelopmental outcomes. This tool predicts binding affinity of chemical compounds to MCT8 using neural network-enhanced molecular docking with Gnina, helping assess potential inhibition risk early in chemical design and drug development.

## Features

- **Web Interface**: User-friendly form for submitting ligands and configuring docking parameters
- **Multiple Input Methods**: SMILES strings, SDF files, or CSV uploads
- **Automated Ligand Preparation**: Hydrogen addition, 3D structure generation, and geometry optimization
- **Neural Network Scoring**: Gnina's CNN-enhanced scoring for improved accuracy
- **Interactive Results**: DataTables with sorting, filtering, and export (CSV, Excel, PDF)
- **3D Molecular Viewer**: Interactive 3D visualization with multi-pose comparison (top 5 poses)
- **Loading Progress Indicator**: Real-time feedback with estimated completion time during docking
- **Risk Assessment**: Automatic categorization of compounds as Likely/Possible/Unlikely inhibitors
- **RESTful API**: Programmatic access for batch processing
- **PDF Reports**: Generate publication-ready docking reports

## Quick Start

### Docker Deployment (Recommended)

The application runs in CPU-only mode using NVIDIA CUDA runtime libraries (no GPU required).

```bash
# Build the Docker image
docker build -t mct8-docking:latest .

# Run the container
docker run -d -p 5000:5000 --name mct8-docking-app \
  -v $(pwd)/results:/usr/src/app/results \
  mct8-docking:latest

# Access the application
open http://localhost:5000

# Or use docker-compose
docker-compose up -d
```

**Container Management:**
```bash
# View logs
docker logs -f mct8-docking-app

# Stop container
docker stop mct8-docking-app

# Restart container
docker start mct8-docking-app

# Remove container
docker rm -f mct8-docking-app
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
   - **CNN Scoring**: none/fast/default (default: none for CPU-only compatibility)

3. **Run Docking**:
   - Click "Run Docking" to execute simulation
   - Or "Generate PDF Report" for downloadable report

4. **Interpret Results**:
   - Binding affinity is reported in kcal/mol (more negative = stronger binding)
   - **Likely Inhibitor**: Binding affinity < -9.0 kcal/mol (High developmental risk)
   - **Possible Inhibitor**: Binding affinity -8.0 to -9.0 kcal/mol (Moderate risk)
   - **Unlikely Inhibitor**: Binding affinity > -8.0 kcal/mol (Low risk)

5. **View in 3D Molecular Viewer**:
   - Click any result row in the table to open the interactive 3D viewer
   - **Mouse Controls**:
     - Left-click + drag: Rotate the view
     - Scroll wheel: Zoom in/out
     - Right-click + drag: Pan the view
   - **Multi-Pose Comparison**: Enable "Show Top 5 Poses" checkbox to overlay the top 5 binding poses in different colors (magenta, cyan, yellow, green, orange)
   - **Toggle Options**:
     - Show/hide protein cartoon backbone
     - Show/hide binding site residues (orange sticks)
   - **Info Panel**: Displays SMILES, binding affinity, CNN score, and risk assessment for the selected pose

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

### Additional API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/receptor` | GET | Returns MCT8 receptor PDB data |
| `/api/binding_site` | GET | Returns binding site PDB data |
| `/api/pose/<id>` | GET | Returns specific docked pose with assessment |
| `/download/sdf` | GET | Downloads complete docking results as SDF file |

#### GET `/api/receptor`

Returns the MCT8 receptor structure in PDB format.

**Response:**
```json
{
  "pdb": "ATOM      1  N   MET A   1...",
  "atoms": 10234,
  "name": "MCT8 Receptor"
}
```

#### GET `/api/binding_site`

Returns the binding site coordinates in PDB format.

**Response:**
```json
{
  "pdb": "HETATM    1  C   SIT E   1...",
  "atoms": 23,
  "name": "Binding Site"
}
```

#### GET `/api/pose/<id>`

Returns a specific docked pose by index (0-based).

**Response:**
```json
{
  "pdb": "ATOM      1  C   LIG A   1...",
  "smiles": "CCO",
  "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
  "affinity": -5.2,
  "cnn_score": 0.65,
  "pose_id": 0,
  "assessment": {
    "category": "Unlikely Inhibitor",
    "color": "#4CAF50",
    "description": "Low inhibition risk (> -8.0 kcal/mol)"
  }
}
```

**cURL Examples:**
```bash
# Get receptor PDB
curl http://localhost:5000/api/receptor

# Get binding site PDB
curl http://localhost:5000/api/binding_site

# Get pose 0 with assessment
curl http://localhost:5000/api/pose/0

# Download all results as SDF
curl -O http://localhost:5000/download/sdf
```

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
6. **Visualization**: 2D structure images, interactive 3D pose viewing

### Scoring

- **Binding Affinity**: Estimated free energy of binding between the ligand and MCT8 receptor (kcal/mol). More negative values indicate stronger binding interactions and higher likelihood of MCT8 inhibition. This is the primary metric for assessing inhibition risk.
- **CNN Score**: Neural network confidence score (0-1) when CNN scoring is enabled
- **Boltzmann Weight**: Thermodynamic probability of each binding pose, calculated from binding affinity and temperature. Higher weights indicate more energetically favorable conformations.

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
- **Python**: 3.10+ (for local development)
- **RAM**: 2GB minimum (4GB recommended for larger batches)
- **Docker**: For containerized deployment (recommended)
- **CPU**: Multi-core recommended (Gnina uses `--cpu 4` flag)

**Note**: The Docker image uses NVIDIA CUDA 12.0 runtime libraries for Gnina compatibility, but **no GPU is required**. The application runs entirely on CPU using Gnina's `--cpu` mode.

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

### CUDA Library Errors (libcudart.so, libcusparse.so, etc.)

The Gnina binary requires CUDA runtime libraries even in CPU-only mode. The Docker image includes these automatically. For local development:

```bash
# The Docker image handles this automatically
# If running locally outside Docker, you may need CUDA runtime libraries
# See: https://developer.nvidia.com/cuda-downloads
```

**Solution**: Use the Docker deployment (recommended) which includes all required CUDA runtime libraries.

### Gnina Download Fails
```bash
# Manually download Gnina v1.3
wget https://github.com/gnina/gnina/releases/download/v1.3/gnina
chmod +x gnina
mkdir -p binaries/
mv gnina binaries/
```

### Docker Build Issues
```bash
# Clear Docker cache
docker system prune -a

# Rebuild without cache
docker build --no-cache -t mct8-docking:latest .
```

### RDKit Import Error
```bash
# Install RDKit conda package (recommended)
conda install -c conda-forge rdkit

# Or use pip (may require system dependencies)
pip install rdkit
```

### Container Won't Start
```bash
# Check logs for errors
docker logs mct8-docking-app

# Verify port 5000 is not in use
lsof -i :5000

# Remove and recreate container
docker rm -f mct8-docking-app
docker run -d -p 5000:5000 --name mct8-docking-app mct8-docking:latest
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
