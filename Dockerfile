FROM python:3.10-slim

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    libxrender1 \
    libexpat1 \
    libxext6 \
    libx11-6 \
    libboost-all-dev \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /usr/src/app

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY . .

# Download Gnina binary
RUN mkdir -p binaries && \
    cd binaries && \
    wget https://github.com/gnina/gnina/releases/download/v1.3/gnina && \
    chmod +x gnina

# Create necessary directories
RUN mkdir -p data temp results

# Expose port
EXPOSE 5000

# Copy and set entrypoint
COPY entrypoint.sh /usr/src/app/entrypoint.sh
RUN chmod +x /usr/src/app/entrypoint.sh

ENTRYPOINT ["/usr/src/app/entrypoint.sh"]
