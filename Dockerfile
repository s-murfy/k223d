# Use an official Python runtime as a parent image
FROM python:3.13-slim

# Set the working directory in the container
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    libglu1-mesa \
    libgl1 \
    libxrender1 \
    libxcursor1 \
    libxinerama1 \
    libxft2 \
    libxrandr2 \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies (copied first so this layer is cached)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy source and build the shared library for the Python API
COPY . /app
RUN make -C /app/code install

# Expose port 8888 for JupyterLab
EXPOSE 8888

# Set environment variable to avoid Jupyter asking for confirmation on launch
ENV JUPYTER_ENABLE_LAB=yes

# Run JupyterLab when the container launches
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--allow-root"]
