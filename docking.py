"""
MCT8 Docking Module
Core docking functionality extracted from Jupyter notebook.
"""

import os
import subprocess
import io
import math
import base64
import logging
from pathlib import Path

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Embedded PDB data
BINDING_SITE_PDB = '''HETATM    1  C   SIT E   1       2.352  -9.145   8.079  1.00 37.89           C
HETATM    2  C   SIT E   1       2.081 -10.116   7.032  1.00 32.58           C
HETATM    3  C   SIT E   1       3.692  -9.231   8.656  1.00 30.81           C
HETATM    4  C   SIT E   1       1.996  -7.771   7.695  1.00 33.02           C
HETATM    5  C   SIT E   1       1.388  -9.481   9.269  1.00 32.44           C
HETATM    6  C   SIT E   1       0.257 -10.187   9.081  1.00 34.12           C
HETATM    7  C   SIT E   1      -0.607 -10.130  10.267  1.00 35.22           C
HETATM    8  C   SIT E   1       0.012 -10.820  11.298  1.00 33.96           C
HETATM    9  C   SIT E   1      -0.849  -8.739  10.746  1.00 31.23           C
HETATM   10  C   SIT E   1      -2.211  -8.596  10.974  1.00 28.57           C
HETATM   11  C   SIT E   1      -0.110  -8.640  12.013  1.00 34.73           C
HETATM   12  C   SIT E   1      -0.740  -7.842  12.929  1.00 28.19           C
HETATM   13  C   SIT E   1      -0.047 -10.055  12.459  1.00 36.34           C
HETATM   14  C   SIT E   1       1.107 -10.324  13.239  1.00 33.75           C
HETATM   15  C   SIT E   1       2.366 -10.297  12.803  1.00 31.31           C
HETATM   16  C   SIT E   1       3.188 -10.640  13.761  1.00 30.74           C
HETATM   17  C   SIT E   1       2.431 -10.891  14.838  1.00 28.96           C
HETATM   18  C   SIT E   1       2.699 -11.294  16.134  1.00 27.96           C
HETATM   19  C   SIT E   1       3.895 -11.526  16.552  1.00 26.80           C
HETATM   20  C   SIT E   1       1.682 -11.452  16.960  1.00 29.89           C
HETATM   21  C   SIT E   1       0.450 -11.222  16.564  1.00 27.00           C
HETATM   22  C   SIT E   1       0.126 -10.841  15.356  1.00 28.72           C
HETATM   23  C   SIT E   1       1.106 -10.681  14.500  1.00 30.95           C'''

# MCT8 receptor PDB loaded from external file
MCT8_RECEPTOR_FILE = "data/mct8_receptor_full.pdb"


def setup_environment():
    """
    Download and setup Gnina (GPU and CPU-capable docking software).
    Gnina can run in CPU-only mode using --cpu flag.

    Returns:
        Path: Path to gnina executable
    """
    gnina_path = Path("binaries/gnina")

    if not gnina_path.exists():
        logger.info("Downloading Gnina v1.3 (CPU-capable)...")
        gnina_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Download Gnina binary for Linux
            subprocess.run([
                "wget",
                "-O", str(gnina_path),
                "https://github.com/gnina/gnina/releases/download/v1.3/gnina"
            ], check=True, capture_output=True)

            gnina_path.chmod(0o755)
            logger.info("Gnina downloaded successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to download Gnina: {e}")
            raise
    else:
        logger.info("Gnina already installed")

    return gnina_path


def create_pdb_files():
    """
    Create PDB files for MCT8 receptor and binding site.

    Returns:
        tuple: Paths to (receptor_pdb, site_pdb)
    """
    data_dir = Path("data")
    data_dir.mkdir(exist_ok=True)

    site_path = data_dir / "binding_site.pdb"
    receptor_path = data_dir / "mct8_receptor.pdb"

    # Write binding site
    with open(site_path, 'w') as f:
        f.write(BINDING_SITE_PDB)

    # Check if full receptor file exists, otherwise create from embedded data
    if Path(MCT8_RECEPTOR_FILE).exists() and receptor_path != Path(MCT8_RECEPTOR_FILE):
        # Copy from full file
        import shutil
        shutil.copy(MCT8_RECEPTOR_FILE, receptor_path)
        logger.info("Using full MCT8 receptor PDB")
    elif not receptor_path.exists():
        # Create placeholder warning - full PDB should be added
        logger.error("Full MCT8 receptor PDB not found!")
        raise FileNotFoundError(f"MCT8 receptor PDB not found at {MCT8_RECEPTOR_FILE}")

    logger.info(f"Created PDB files: {site_path}, {receptor_path}")
    return receptor_path, site_path


def validate_smiles(smiles_list):
    """
    Validate SMILES strings using RDKit.

    Args:
        smiles_list (list): List of SMILES strings

    Returns:
        dict: {'valid': [...], 'invalid': [...]}
    """
    valid = []
    invalid = []

    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles.strip())
            if mol is not None:
                valid.append(smiles.strip())
            else:
                invalid.append(smiles.strip())
        except Exception as e:
            logger.warning(f"Error validating SMILES '{smiles}': {e}")
            invalid.append(smiles.strip())

    return {'valid': valid, 'invalid': invalid}


def prepare_ligands(smiles_list, output_file="ligands.sdf",
                   add_hydrogens=True, generate_3d=True, optimize=True):
    """
    Prepare ligands from SMILES strings.

    Args:
        smiles_list (list): List of SMILES strings
        output_file (str): Output SDF filename
        add_hydrogens (bool): Add hydrogen atoms
        generate_3d (bool): Generate 3D coordinates
        optimize (bool): Optimize structure with MMFF

    Returns:
        str: Path to output SDF file
    """
    sdf_writer = Chem.SDWriter(output_file)

    for i, smiles in enumerate(smiles_list):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(f"Skipping invalid SMILES: {smiles}")
                continue

            # Add hydrogens
            if add_hydrogens:
                try:
                    mol = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
                except:
                    logger.warning(f"Mol({i}): Add hydrogens failed!")

            # Generate 3D structure
            if generate_3d:
                try:
                    embed_params = AllChem.ETKDGv3()
                    embed_params.randomSeed = 0xf00d
                    AllChem.EmbedMolecule(mol, embed_params)
                except:
                    logger.warning(f"Mol({i}): 3D structure generation failed!")

            # Optimize structure
            if optimize:
                try:
                    AllChem.MMFFOptimizeMolecule(mol)
                except:
                    logger.warning(f"Mol({i}): Structure optimization failed!")

            # Write to SDF
            try:
                sdf_writer.write(mol)
            except:
                logger.warning(f"Mol({i}): Structure save failed!")

        except Exception as e:
            logger.error(f"Error processing SMILES '{smiles}': {e}")

    sdf_writer.close()
    logger.info(f"Prepared ligands saved to {output_file}")
    return output_file


def run_docking(ligand_file, receptor_file, site_file,
                num_modes=3, exhaustiveness=8, autobox_add=6.0,
                cnn="none", output_file="output.sdf", log_file="output.txt"):
    """
    Run Gnina molecular docking with CPU-only mode support.

    Args:
        ligand_file (str): Path to ligand SDF file
        receptor_file (str): Path to receptor PDB file
        site_file (str): Path to binding site PDB file (for autobox)
        num_modes (int): Number of binding modes to generate
        exhaustiveness (int): Search exhaustiveness
        autobox_add (float): Autobox expansion in Angstroms
        cnn (str): CNN scoring model ('none', 'fast', 'default')
        output_file (str): Output file
        log_file (str): Log file

    Returns:
        dict: {'output_file': str, 'log_file': str, 'success': bool, 'error': str or None}
    """
    gnina_path = Path("binaries/gnina")

    if not gnina_path.exists():
        return {
            'success': False,
            'error': "Gnina not found. Run setup_environment() first."
        }

    # Build Gnina command
    cmd = [
        str(gnina_path),
        "-r", receptor_file,
        "-l", ligand_file,
        "--autobox_ligand", site_file,
        "--autobox_add", str(autobox_add),
        "-o", output_file,
        "--num_modes", str(num_modes),
        "--exhaustiveness", str(exhaustiveness),
        "--cpu", "4"  # CPU-only mode (no CUDA required)
    ]

    # Handle CNN scoring
    if cnn == "none":
        cmd.extend(["--cnn_scoring", "none"])
    elif cnn == "fast":
        cmd.extend(["--cnn", "crossdock_default2018"])
    elif cnn == "default":
        cmd.extend(["--cnn", "default2017"])

    # Add log file
    if log_file:
        cmd.extend(["--log", log_file])

    logger.info(f"Running Gnina docking (CPU mode): {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode == 0:
            logger.info("Docking completed successfully")
            return {
                'success': True,
                'output_file': output_file,
                'log_file': log_file,
                'error': None,
                'stderr': result.stderr
            }
        else:
            logger.error(f"Docking failed: {result.stderr}")
            return {
                'success': False,
                'output_file': None,
                'log_file': None,
                'error': result.stderr
            }

    except subprocess.TimeoutExpired:
        logger.error("Docking timeout (10 minutes)")
        return {
            'success': False,
            'error': "Docking timeout after 10 minutes"
        }
    except Exception as e:
        logger.error(f"Docking error: {e}")
        return {
            'success': False,
            'error': str(e)
        }


def process_results(sdf_file):
    """
    Process docking results from SDF file.

    Args:
        sdf_file (str): Path to output SDF file

    Returns:
        pd.DataFrame: Results with affinity, CNN scores, Boltzmann weights
    """
    try:
        poses = Chem.SDMolSupplier(sdf_file, removeHs=False)
        data = []

        for p_idx, pose in enumerate(poses):
            if pose is None:
                continue

            props = pose.GetPropsAsDict()
            result = {
                'pose_index': p_idx,
                'smiles': Chem.MolToSmiles(pose),
                'inchikey': Chem.inchi.MolToInchiKey(pose),
                'affinity': props.get('minimizedAffinity', None),
                'cnn_score': props.get('CNNscore', None),
                'cnn_affinity': props.get('CNNaffinity', None)
            }

            # Calculate Boltzmann partition function component
            if result['affinity'] is not None:
                result['q_part'] = math.exp(-1 * result['affinity'] / 1.987204259 * 0.3)
            else:
                result['q_part'] = 0

            data.append(result)

        df = pd.DataFrame(data)

        # Calculate Boltzmann weights
        if len(df) > 0 and 'q_part' in df.columns:
            for inchikey in df['inchikey'].unique():
                mask = df['inchikey'] == inchikey
                total_q = df.loc[mask, 'q_part'].sum()
                if total_q > 0:
                    df.loc[mask, 'boltzmann_weight'] = df.loc[mask, 'q_part'] / total_q

        logger.info(f"Processed {len(df)} docking poses")
        return df

    except Exception as e:
        logger.error(f"Error processing results: {e}")
        return pd.DataFrame()


def smiles_to_image(smiles, img_size=(300, 300)):
    """
    Convert SMILES to PNG image (base64 data URI).

    Args:
        smiles (str): SMILES string
        img_size (tuple): Image size (width, height)

    Returns:
        str: Base64 data URI or None
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        img = Draw.MolToImage(mol, size=img_size)

        # Convert to base64 data URI
        buffer = io.BytesIO()
        img.save(buffer, format='PNG')
        img_str = base64.b64encode(buffer.getvalue()).decode()

        return f"data:image/png;base64,{img_str}"

    except Exception as e:
        logger.warning(f"Error converting SMILES to image: {e}")
        return None


def assess_inhibition(affinity):
    """
    Assess inhibition risk based on affinity score.

    Args:
        affinity (float): Binding affinity in kcal/mol

    Returns:
        dict: {'category': str, 'color': str, 'description': str}
    """
    if affinity is None:
        return {
            'category': 'Unknown',
            'color': '#999999',
            'description': 'No affinity data'
        }

    if affinity < -9.0:
        return {
            'category': 'Strong Inhibitor',
            'color': '#E6007E',  # Magenta
            'description': 'High inhibition risk (< -9.0 kcal/mol)'
        }
    elif affinity < -8.0:
        return {
            'category': 'Moderate Inhibitor',
            'color': '#FF9500',  # Orange
            'description': 'Moderate inhibition risk (-8.0 to -9.0 kcal/mol)'
        }
    else:
        return {
            'category': 'Weak/Non-Inhibitor',
            'color': '#4CAF50',  # Green
            'description': 'Low inhibition risk (> -8.0 kcal/mol)'
        }


if __name__ == "__main__":
    # Test module functionality
    print("MCT8 Docking Module")
    print("=" * 50)

    # Setup
    print("\n1. Setting up environment...")
    setup_environment()

    # Create PDB files
    print("\n2. Creating PDB files...")
    receptor, site = create_pdb_files()

    # Test SMILES validation
    print("\n3. Testing SMILES validation...")
    test_smiles = [
        "CCO",  # Ethanol - valid
        "INVALID",  # Invalid
        "c1ccccc1"  # Benzene - valid
    ]
    validation = validate_smiles(test_smiles)
    print(f"Valid: {validation['valid']}")
    print(f"Invalid: {validation['invalid']}")

    print("\nModule test complete!")
