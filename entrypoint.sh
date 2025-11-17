#!/bin/sh
# MCT8 Docking Docker Entrypoint
echo "Starting MCT8 Docking Flask Application..."
echo "CPU-only mode enabled (--cpu 4)"
echo "Server will be available at http://0.0.0.0:5000"
exec python app.py
