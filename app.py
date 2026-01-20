"""
MCT8 Docking Flask Application
Web interface for molecular docking simulations.
"""

import os
import io
import logging
from pathlib import Path
from flask import Flask, render_template, request, jsonify, send_file, flash, redirect
from flask_cors import CORS

import pandas as pd
from werkzeug.utils import secure_filename
from rdkit import Chem

import docking

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize Flask app
app = Flask(__name__)
app.secret_key = os.environ.get('SECRET_KEY', 'dev-secret-key-change-in-production')
CORS(app)

# Configuration
ALLOWED_EXTENSIONS = {'sdf', 'smi', 'csv', 'txt'}
MAX_SMILES_INPUT = 100


def allowed_file(filename):
    """Check if file extension is allowed."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def validate_docking_params(params):
    """
    Validate docking parameters.

    Args:
        params (dict): Dictionary with num_modes, exhaustiveness, autobox_add, cnn

    Returns:
        list: List of error messages (empty if valid)
    """
    errors = []

    num_modes = params.get('num_modes')
    if num_modes is not None:
        if not isinstance(num_modes, int) or not 1 <= num_modes <= 50:
            errors.append("num_modes must be an integer between 1 and 50")

    exhaustiveness = params.get('exhaustiveness')
    if exhaustiveness is not None:
        if not isinstance(exhaustiveness, int) or not 2 <= exhaustiveness <= 128:
            errors.append("exhaustiveness must be an integer between 2 and 128")

    autobox_add = params.get('autobox_add')
    if autobox_add is not None:
        if not isinstance(autobox_add, (int, float)) or not 0 <= autobox_add <= 10:
            errors.append("autobox_add must be a number between 0 and 10")

    cnn = params.get('cnn')
    if cnn is not None:
        if cnn not in ['none', 'fast', 'default']:
            errors.append("cnn must be 'none', 'fast', or 'default'")

    return errors


@app.route('/')
def home():
    """Render home page with docking form."""
    return render_template('index.html')


@app.route('/dock', methods=['POST'])
def dock():
    """
    Handle docking request from web form.

    Form data:
        - smiles: Comma-separated SMILES strings
        - file: SDF/SMI file upload
        - add_h: Add hydrogens (checkbox)
        - gen_3d: Generate 3D structure (checkbox)
        - optimize: Optimize structure (checkbox)
        - num_modes: Number of binding modes (1-50)
        - exhaustiveness: Search exhaustiveness (2-128)
        - autobox_add: Search box expansion (0-10)
        - cnn: CNN scoring model ('none', 'fast', 'default')
        - action: 'dock' or 'report'

    Returns:
        Rendered template with results table
    """
    try:
        # Initialize environment
        logger.info("Setting up docking environment...")
        docking.setup_environment()
        receptor, site = docking.create_pdb_files()

        # Get SMILES from text input or file
        smiles_list = []

        # Text input
        smiles_input = request.form.get('smiles', '').strip()
        if smiles_input:
            smiles_list.extend([s.strip() for s in smiles_input.split(',') if s.strip()])

        # File upload
        if 'file' in request.files:
            file = request.files['file']
            if file and file.filename and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                ext = filename.rsplit('.', 1)[1].lower()

                if ext == 'sdf':
                    # Save temporarily and read
                    temp_path = f"temp_{filename}"
                    file.save(temp_path)
                    # Will be processed directly by prepare_ligands
                    ligand_file = temp_path
                    logger.info(f"Uploaded SDF file: {filename}")
                elif ext in ['smi', 'csv', 'txt']:
                    # Read SMILES from file
                    content = file.read().decode('utf-8')
                    if ext == 'csv':
                        # Try to read as CSV with 'SMILES' column
                        df = pd.read_csv(io.StringIO(content))
                        if 'SMILES' in df.columns:
                            smiles_list.extend(df['SMILES'].dropna().tolist())
                        else:
                            flash('CSV file must have a "SMILES" column', 'error')
                            return render_template('index.html', error='CSV file must have a "SMILES" column')
                    else:
                        # Read line by line
                        smiles_list.extend([s.strip() for s in content.split('\n') if s.strip()])

        if not smiles_list and not ('file' in request.files and file and file.filename.endswith('.sdf')):
            flash('Please provide SMILES strings or upload a file', 'error')
            return render_template('index.html', error='No input provided')

        # Limit number of SMILES
        if len(smiles_list) > MAX_SMILES_INPUT:
            flash(f'Maximum {MAX_SMILES_INPUT} SMILES allowed. First {MAX_SMILES_INPUT} will be processed.', 'warning')
            smiles_list = smiles_list[:MAX_SMILES_INPUT]

        # Validate SMILES (only if we have SMILES, not SDF)
        if smiles_list:
            logger.info(f"Validating {len(smiles_list)} SMILES...")
            validation = docking.validate_smiles(smiles_list)
            valid_smiles = validation['valid']
            invalid_smiles = validation['invalid']

            if invalid_smiles:
                flash(f'Invalid SMILES ({len(invalid_smiles)}): {", ".join(invalid_smiles[:5])}...', 'warning')

            if not valid_smiles:
                return render_template('index.html', error='No valid SMILES provided',
                                     invalid_smiles=invalid_smiles)

            # Prepare ligands
            logger.info(f"Preparing {len(valid_smiles)} ligands...")
            ligand_file = docking.prepare_ligands(
                valid_smiles,
                output_file="ligands.sdf",
                add_hydrogens=request.form.get('add_h') == 'on',
                generate_3d=request.form.get('gen_3d') == 'on',
                optimize=request.form.get('optimize') == 'on'
            )

        # Get docking parameters
        params = {
            'num_modes': int(request.form.get('num_modes', 3)),
            'exhaustiveness': int(request.form.get('exhaustiveness', 8)),
            'autobox_add': float(request.form.get('autobox_add', 6.0)),
            'cnn': request.form.get('cnn', 'none')  # Default to 'none' (CPU-only)
        }

        # Validate parameters
        validation_errors = validate_docking_params(params)
        if validation_errors:
            return render_template('index.html', error=f"Invalid parameters: {'; '.join(validation_errors)}")

        logger.info(f"Running docking with parameters: {params}")

        # Run docking
        result = docking.run_docking(
            ligand_file=ligand_file,
            receptor_file=str(receptor),
            site_file=str(site),
            **params
        )

        if not result['success']:
            return render_template('index.html', error=f"Docking failed: {result['error']}")

        # Process results
        logger.info("Processing docking results...")
        df = docking.process_results(result['output_file'])

        if df.empty:
            return render_template('index.html', error='No docking results generated')

        # Add molecular images and assessments
        df['structure_img'] = df['smiles'].apply(docking.smiles_to_image)
        df['assessment'] = df['affinity'].apply(lambda x: docking.assess_inhibition(x) if pd.notna(x) else {})

        # Check if user requested PDF report
        if request.form.get('action') == 'report':
            # Generate PDF report
            logger.info("Generating PDF report...")
            pdf_buffer = generate_pdf_report(df, params)
            return send_file(
                pdf_buffer,
                mimetype='application/pdf',
                as_attachment=True,
                download_name='mct8_docking_report.pdf'
            )

        # Render results
        results_html = df.to_dict('records')
        return render_template('index.html',
                             results=results_html,
                             params=params,
                             num_poses=len(df))

    except Exception as e:
        logger.error(f"Error in docking: {e}", exc_info=True)
        return render_template('index.html', error=f'An error occurred: {str(e)}')


@app.route('/api', methods=['POST'])
def api_dock():
    """
    RESTful API endpoint for docking.

    JSON Request:
        {
            "smiles": ["SMILES1", "SMILES2", ...],
            "params": {
                "num_modes": 3,
                "exhaustiveness": 8,
                "autobox_add": 6.0,
                "cnn": "fast"
            },
            "format": "json|csv|sdf"
        }

    Returns:
        JSON/CSV/SDF with docking results
    """
    try:
        data = request.get_json()

        if not data or 'smiles' not in data:
            return jsonify({'error': 'Missing SMILES data'}), 400

        smiles_list = data['smiles']
        if not isinstance(smiles_list, list):
            return jsonify({'error': 'SMILES must be a list'}), 400

        if len(smiles_list) > MAX_SMILES_INPUT:
            return jsonify({'error': f'Maximum {MAX_SMILES_INPUT} SMILES allowed'}), 400

        # Validate SMILES
        validation = docking.validate_smiles(smiles_list)
        valid_smiles = validation['valid']
        invalid_smiles = validation['invalid']

        if not valid_smiles:
            return jsonify({
                'error': 'No valid SMILES provided',
                'invalid_smiles': invalid_smiles
            }), 400

        # Initialize
        docking.setup_environment()
        receptor, site = docking.create_pdb_files()

        # Prepare ligands
        ligand_file = docking.prepare_ligands(valid_smiles)

        # Get parameters
        params = data.get('params', {})
        docking_params = {
            'num_modes': params.get('num_modes', 3),
            'exhaustiveness': params.get('exhaustiveness', 8),
            'autobox_add': params.get('autobox_add', 6.0),
            'cnn': params.get('cnn', 'none')  # Default to 'none' (CPU-only, no CUDA required)
        }

        # Validate parameters
        validation_errors = validate_docking_params(docking_params)
        if validation_errors:
            return jsonify({'error': f"Invalid parameters: {'; '.join(validation_errors)}"}), 400

        # Run docking
        result = docking.run_docking(
            ligand_file=ligand_file,
            receptor_file=str(receptor),
            site_file=str(site),
            **docking_params
        )

        if not result['success']:
            return jsonify({'error': result['error']}), 500

        # Process results
        df = docking.process_results(result['output_file'])

        if df.empty:
            return jsonify({'error': 'No results generated'}), 500

        # Add assessments
        df['assessment'] = df['affinity'].apply(
            lambda x: docking.assess_inhibition(x)['category'] if pd.notna(x) else 'Unknown'
        )

        # Format output
        output_format = data.get('format', 'json')

        if output_format == 'csv':
            csv_buffer = io.StringIO()
            df.to_csv(csv_buffer, index=False)
            return csv_buffer.getvalue(), 200, {'Content-Type': 'text/csv'}

        elif output_format == 'sdf':
            return send_file(
                result['output_file'],
                mimetype='chemical/x-mdl-sdfile',
                as_attachment=True,
                download_name='docking_results.sdf'
            )

        else:  # json
            results = df.to_dict('records')
            return jsonify({
                'success': True,
                'num_poses': len(results),
                'invalid_smiles': invalid_smiles,
                'results': results
            })

    except Exception as e:
        logger.error(f"API error: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500


@app.route('/api/receptor')
def api_receptor():
    """Get MCT8 receptor PDB file."""
    try:
        receptor_path = Path("data/mct8_receptor.pdb")
        if not receptor_path.exists():
            return jsonify({'error': 'Receptor PDB not found'}), 404

        with open(receptor_path, 'r') as f:
            pdb_data = f.read()

        atom_count = len([line for line in pdb_data.split('\n') if line.startswith('ATOM')])

        return jsonify({
            'pdb': pdb_data,
            'atoms': atom_count,
            'name': 'MCT8 Receptor'
        })
    except Exception as e:
        logger.error(f"Error loading receptor: {e}")
        return jsonify({'error': str(e)}), 500


@app.route('/api/binding_site')
def api_binding_site():
    """Get binding site PDB file."""
    try:
        site_path = Path("data/binding_site.pdb")
        if not site_path.exists():
            return jsonify({'error': 'Binding site PDB not found'}), 404

        with open(site_path, 'r') as f:
            pdb_data = f.read()

        atom_count = len([line for line in pdb_data.split('\n') if line.startswith('HETATM')])

        return jsonify({
            'pdb': pdb_data,
            'atoms': atom_count,
            'name': 'Binding Site'
        })
    except Exception as e:
        logger.error(f"Error loading binding site: {e}")
        return jsonify({'error': str(e)}), 500


@app.route('/api/pose/<int:pose_id>')
def api_pose(pose_id):
    """Get PDB data for a specific docked pose."""
    try:
        output_file = "output.sdf"
        if not Path(output_file).exists():
            return jsonify({'error': 'No docking results found. Run docking first.'}), 404

        # Read SDF file
        supplier = Chem.SDMolSupplier(output_file, removeHs=False)
        poses = list(supplier)

        # Validate pose_id
        if pose_id < 0 or pose_id >= len(poses):
            return jsonify({'error': f'Invalid pose_id. Must be 0-{len(poses)-1}'}), 400

        pose = poses[pose_id]
        if pose is None:
            return jsonify({'error': f'Pose {pose_id} could not be loaded'}), 500

        # Convert to PDB format
        pdb_block = Chem.MolToPDBBlock(pose)
        pdb_block = pdb_block.replace("UNL", "LIG")  # Replace generic residue name

        # Get properties
        props = pose.GetPropsAsDict()
        affinity = props.get('minimizedAffinity', None)
        cnn_score = props.get('CNNscore', None)

        # Generate response
        result = {
            'pdb': pdb_block,
            'smiles': Chem.MolToSmiles(pose),
            'inchikey': Chem.inchi.MolToInchiKey(pose),
            'affinity': affinity,
            'cnn_score': cnn_score,
            'pose_id': pose_id
        }

        # Add assessment
        if affinity is not None:
            result['assessment'] = docking.assess_inhibition(affinity)

        return jsonify(result)

    except Exception as e:
        logger.error(f"Error loading pose {pose_id}: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500


@app.route('/download/sdf')
def download_sdf():
    """Download complete output.sdf file with all docked poses."""
    try:
        output_file = Path("output.sdf")
        if not output_file.exists():
            flash('No docking results found. Please run docking first.', 'error')
            return redirect('/')

        return send_file(
            output_file,
            mimetype='chemical/x-mdl-sdfile',
            as_attachment=True,
            download_name='mct8_docking_results.sdf'
        )
    except Exception as e:
        logger.error(f"Error downloading SDF: {e}")
        flash(f'Error downloading SDF file: {str(e)}', 'error')
        return redirect('/')


def generate_pdf_report(df, params):
    """
    Generate PDF report of docking results.

    Args:
        df: DataFrame with results
        params: Docking parameters dict

    Returns:
        BytesIO buffer with PDF
    """
    from reportlab.lib.pagesizes import letter, A4
    from reportlab.lib import colors
    from reportlab.lib.units import inch
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
    from reportlab.lib.styles import getSampleStyleSheet

    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    elements = []
    styles = getSampleStyleSheet()

    # Title
    title = Paragraph("<b>MCT8 Docking Results Report</b>", styles['Title'])
    elements.append(title)
    elements.append(Spacer(1, 0.2 * inch))

    # Parameters
    params_text = f"""<b>Docking Parameters:</b><br/>
    Number of Modes: {params.get('num_modes', 'N/A')}<br/>
    Exhaustiveness: {params.get('exhaustiveness', 'N/A')}<br/>
    Search Box Expansion: {params.get('autobox_add', 'N/A')} Å<br/>
    CNN Scoring: {params.get('cnn', 'N/A')}
    """
    elements.append(Paragraph(params_text, styles['Normal']))
    elements.append(Spacer(1, 0.3 * inch))

    # Results table
    table_data = [['SMILES', 'Affinity (kcal/mol)', 'CNN Score', 'Assessment']]
    for _, row in df.head(20).iterrows():  # Limit to first 20 results
        assessment = docking.assess_inhibition(row.get('affinity'))
        table_data.append([
            row.get('smiles', '')[:40] + '...' if len(row.get('smiles', '')) > 40 else row.get('smiles', ''),
            f"{row.get('affinity', 'N/A'):.2f}" if pd.notna(row.get('affinity')) else 'N/A',
            f"{row.get('cnn_score', 'N/A'):.3f}" if pd.notna(row.get('cnn_score')) else 'N/A',
            assessment['category']
        ])

    table = Table(table_data)
    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black)
    ]))

    elements.append(table)

    # Build PDF
    doc.build(elements)
    buffer.seek(0)
    return buffer


if __name__ == '__main__':
    debug_mode = os.environ.get('FLASK_DEBUG', 'false').lower() == 'true'
    app.run(host='0.0.0.0', port=5000, debug=debug_mode)
